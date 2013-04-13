import json
import re
import sqlite3
from chilin2.helpers import JinjaTemplateCommand, template_dump, r_exec, json_dump, json_load

def _python_fastqc_parse(input, output=None, param=None):
    data = open(input).readlines()
    sequence_length = 0
    quality_dict = {}
    in_seq_quality_section = False
    for line in data:
        if re.search(r"^Sequence length", line):
            assert sequence_length == 0
            sequence_length = int(re.findall(r"^Sequence length\t(\d+)", line)[0])
        elif re.search(r"^>>Per sequence quality", line):
            assert not in_seq_quality_section
            in_seq_quality_section = True
            continue

        if re.search(r"^>>END_MODULE", line) and in_seq_quality_section:
            in_seq_quality_section = False

        if (not line.startswith("#")) and in_seq_quality_section:
            sequence_quality = re.findall("^(\w+)\t(\w+)", line)[0]
            quality_dict[sequence_quality[0]] = float(sequence_quality[1])
    total = sum(quality_dict.values())
    n = 0
    for item in sorted(quality_dict.items(), key=lambda e: e[0], reverse=True):
        n = n + item[1]
        if n / total > 0.5:
            median = int(item[0])
            break
    return {"sequence_length": sequence_length,
            "median": median}


def stat_fastqc(input={"db": "", "fastqc_summaries": [], "template": ""},
                output={"R": "", "json": "", "pdf": ""},
                param={"ids": [], "id": ""}):
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    stat = json_dict["stat"]

    quality_medians = []

    for a_summary, a_id in zip(input["fastqc_summaries"], param["ids"]):
        parsed = _python_fastqc_parse(input=a_summary)

        stat[a_id] = {}
        stat[a_id]["median"] = parsed["median"]
        stat[a_id]["cutoff"] = 25
        stat[a_id]['judge'] = "Pass" if parsed["median"] > 25 else "Fail"
        stat[a_id]["sequence_length"] = parsed["sequence_length"]

        quality_medians.append(parsed["median"])

    # The table of fastqc_summary that will be used for rendering
    # Col 1: sample ID
    # Col 2: sequence length
    # Col 3: median of sequence quality

    qc_db = sqlite3.connect(input["db"]).cursor()
    qc_db.execute("SELECT median_quality FROM fastqc_info")
    history_data = [float(i[0]) for i in qc_db.fetchall()]

    fastqc_dist_r = JinjaTemplateCommand(
        template=input["template"],
        param={'historic_data': history_data,
               'current_data': quality_medians,
               'ids': param["ids"],
               'cutoff': 25,
               'main': 'Sequence Quality Score Cumulative Percentage',
               'xlab': 'sequence quality score',
               'ylab': 'fn(sequence quality score)',
               "need_smooth_curve": True,

               "pdf": output["pdf"],
               "render_dump": output["R"]})

    template_dump(fastqc_dist_r)
    r_exec(fastqc_dist_r)

    json_dump(json_dict)



def latex_fastqc(input, output, param):
    json_dict = json_load(input["json"])

    fastqc_summary = []
    stat = json_dict["stat"]
    for sample in stat:
        fastqc_summary.append([sample, stat[sample]["sequence_length"], stat[sample]["median"]])

    latex = JinjaTemplateCommand(
        template=input["template"],
        param={"section_name": "sequence_quality",
               "path": json_dict["output"]["pdf"],
               "fastqc_table": fastqc_summary,
               "fastqc_graph": json_dict["output"]["pdf"],
               'prefix_dataset_id': stat.keys(),

               "render_dump": output["latex"]})
    template_dump(latex)
