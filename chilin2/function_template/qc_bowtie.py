import json
import re
import sqlite3
from chilin2.helpers import JinjaTemplateCommand, template_dump, r_exec, json_load


def underline_to_space(x):
    if type(x) == str:
        return x.replace("_", " ")
    return x


def _bowtie_summary_parse(input=""):
    summary_content = open(input).read()
    print(summary_content)
    print("_" * 100)
    total_reads = int(re.findall(r"# reads processed: (.*)\n", summary_content)[0])

    # WARN: this is `unique mappable reads` as we set `-m 1`
    mappable_reads = int(re.findall(r"# reads with at least one reported alignment: (.*) \(", summary_content)[0])

    mappable_rate = float(mappable_reads) / float(total_reads)

    return {"total_reads": total_reads, "mappable_reads": mappable_reads, "mappable_rate": mappable_rate}


def stat_bowtie(input={"bowtie_summaries": [], "db":"", "template": ""},
                output={"json": "", "R": "", "pdf": ""},
                param={"sams": [], }):
    """ sams = [{'name1':a, 'total1': 5...}, {'name2':c, 'total2': 3...}...] **args """

    # unique location is in text_macs2_summary part
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}

    db = sqlite3.connect(input["db"]).cursor()
    db.execute("select map_ratio from mapping")
    historyData = [str(i[0]) for i in (db.fetchall())]
    bowtie_summaries = {"total_reads": [],
                        "mappable_reads": [],
                        "mappable_rate": []}

    for summary, sam in zip(input["bowtie_summaries"], param["sams"]):
        json_dict["stat"][sam] = _bowtie_summary_parse(summary)

    mappable_rates = [json_dict["stat"][i]["mappable_rate"] for i in json_dict["stat"]]

    mappable_rate_R = JinjaTemplateCommand(template=input["template"],
        param={'historic_data': historyData,
               'current_data': mappable_rates,
               'ids': param["sams"],
               'cutoff': 0.5,
               'main': 'Unique mapped rate',
               'xlab': 'Unique mapped rate',
               'ylab': 'fn(Unique mapped rate)',
               "need_smooth_curve": True,

               "render_dump": output["R"],
               "pdf": output["pdf"], })
    template_dump(mappable_rate_R)
    r_exec(mappable_rate_R)

    with open(output["json"], "w") as f:
        json.dump(json_dict, f, indent=4)


def latex_bowtie(input, output, param):
    json_dict = json_load(input["json"])

    basic_map_table = []
    for sam in json_dict["stat"]:
        basic_map_table.append([underline_to_space(sam),
                                json_dict["stat"][sam]["total_reads"],
                                json_dict["stat"][sam]["mappable_reads"],
                                json_dict["stat"][sam]["mappable_rate"]])

    latex = JinjaTemplateCommand(
        name="mapping quality",
        template=input["template"],
        param={"section_name": "bowtie",
               "basic_map_table": basic_map_table,
               "mappable_ratio_graph": json_dict["output"]["pdf"],

               "render_dump": output["latex"]})

    template_dump(latex)
