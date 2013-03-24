import re
import sqlite3
from pyflow.command import ThrowableShellCommand
from chilin2.jinja_template_render import JinjaTemplateCommand, write_and_run_Rscript, write_into


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

    return {"total_reads": total_reads, "mappable_reads": mappable_reads, "mappable_rate": mappable_rate}


def qc_bowtie_summary_draw(input={"all_bowtie_summary": "", "db": "", "R_template": "", "latex_template": ""},
                           output={"rfile": "", "latex_section": "", "pdf": ""},
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

    write_and_run_Rscript(mappable_rate_R, output["rfile"])

    # basic mappable table, two layer list
    # COL 1 total reads
    # COL 2 mappable reads
    # COL 3 mappable rates
    basic_map_table = []
    for l in range(len(current_data_bowtie_summary)):
        basic_map_table.append(
            [ current_data_bowtie_summary[i][l] for i in current_data_bowtie_summary ])
        
    mapping_quality_latex = JinjaTemplateCommand(
        name="mapping quality",
        template=input["latex_template"],
        param={"section_name": "bowtie",
               "basic_map_table": basic_map_table,
               "mappable_ratio_graph": output["pdf"],
               })
    print(output["latex_section"])
    
    write_into(mapping_quality_latex, output['latex_section'])
    return {}
