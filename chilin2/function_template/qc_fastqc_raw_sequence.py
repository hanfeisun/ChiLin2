import re
import sqlite3
from samflow.command import ThrowableShellCommand
from chilin2.jinja_template_render import JinjaTemplateCommand, write_and_run_Rscript, write_into

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


def python_fastqc_dist_draw(input={"db": "", "fastqc_summary_list": [], "R_template": "", "latex_template": ""},
                            output={"rfile": "", "latex_section": "", "pdf": ""},
                            param={"ids": [], "id": ""}):
    seq_lengths = []
    quality_medians = []


    for a_summary in input["fastqc_summary_list"]:
        parsed_result = _python_fastqc_parse(input=a_summary)
        quality_medians.append(parsed_result["median"])

    # The table of fastqc_summary that will be used for rendering
    # Col 1: sample ID
    # Col 2: sequence length
    # Col 3: median of sequence quality
    fastqc_summary = []
    for i in zip(param["ids"], seq_lengths, quality_medians):
        fastqc_summary.append(list(i))

    qc_db = sqlite3.connect(input["db"]).cursor()
    qc_db.execute("SELECT median_quality FROM fastqc_info")
    history_data = [float(i[0]) for i in qc_db.fetchall()]

    fastqc_dist_R = JinjaTemplateCommand(
        template=input["R_template"],
        param={'historic_data': history_data,
               'current_data': quality_medians,
               'ids': param["ids"],
               'cutoff': 25,
               'main': 'Sequence Quality Score Cumulative Percentage',
               'xlab': 'sequence quality score',
               'ylab': 'fn(sequence quality score)',
               "pdf": output["pdf"],
               "need_smooth_curve": True})

    write_and_run_Rscript(fastqc_dist_R, output["rfile"])
    print(fastqc_summary)
    sequence_quality_latex = JinjaTemplateCommand(
        name = "sequence quality",
        template=input['latex_template'],
        param={"section_name": "sequence_quality",
               "path": output["pdf"],
               "fastqc_table": fastqc_summary,
               "fastqc_graph": output["pdf"],
               'prefix_dataset_id': param['id']})
    write_into(sequence_quality_latex, output["latex_section"])

    return {}


