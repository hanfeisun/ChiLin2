import os
import re
import sqlite3
from pyflow.command import ThrowableShellCommand
from pyflow.helper import fetch

# TODO (Hanfei): group here doesn't work now
#digit_format = lambda num,sep=",": sep[::-1].join(group(str(num)[::-1],3))[::-1]
from chilin2.jinja_template_command import JinjaTemplateCommand

def underline_to_space(x):
    if type(x) == str:
        return x.replace("_", " ")
    return x


def _qc_bowtie_summary_parse(input=""):
    summary_content = open(input).read()
    print(summary_content)
    print("_"*100)
    total_reads = int(re.findall(r"# reads processed: (.*)\n", summary_content)[0])
    
    # WARN: this is `unique mappable reads` as we set `-m 1`
    mappable_reads = int(re.findall(r"# reads with at least one reported alignment: (.*) \(", summary_content)[0])

    mappable_rate = float(mappable_reads) / float(total_reads)

    return {"total_reads": total_reads,
            "mappable_reads": mappable_reads,
            "mappable_rate": mappable_rate}


def qc_bowtie_summary_draw(input={"all_bowtie_summary": "", "db": "", "R_template": "", "Latex_template": ""},
                           output={"rfile": "", "pdf": ""},
                           param={"ids": []}):
    db = sqlite3.connect(input["db"]).cursor()
    db.execute("select map_ratio from mapping")
    historyData = [str(i[0]) for i in (db.fetchall())]
    current_data_bowtie_summary = {"total_reads": [],
                                   "mappable_reads": [],
                                   "mappable_rate": []}

    for a_summary in input["all_bowtie_summary"]:
        parsed_summary = _qc_bowtie_summary_parse(a_summary)
        for k, v in parsed_summary.items():
            current_data_bowtie_summary[k].append(v)

    mappable_rate_R = JinjaTemplateCommand(template=input["R_template"],
        param={'historic_data': historyData,
               'current_data': current_data_bowtie_summary["mappable_rate"],
               'ids': param["ids"],
               'cutoff': 0.5,
               'main': 'Unique mapped rate',
               'xlab': 'Unique mapped rate',
               'ylab': 'fn(Unique mapped rate)',
               "pdf": output["pdf"],
               "need_smooth_curve": True})

    mappable_rate_R.invoke()
    with open(output["rfile"], 'w') as f:
        f.write(mappable_rate_R.result)

    ThrowableShellCommand(template = 'Rscript {input}', input= output["rfile"]).invoke()
    return {}
