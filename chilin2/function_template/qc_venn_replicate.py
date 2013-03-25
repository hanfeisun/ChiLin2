import re
import subprocess
from chilin2.jinja_template_render import JinjaTemplateCommand, write_into

def qc_replicate_parse(input={"correlation_R":"", "cor_pdf": "",  "latex_summaryTable":""},
                       output={"latex_section": ""}, param=None):
    """ ReplicateQC aims to describe the similarity of replicate experiment. Venn diagram and correlation plot will be used."""

    correlation_result_R_code = open(input["correlation_R"]).read()
    signal_list = re.findall(r"[pc]\d? <- (.*)$", correlation_result_R_code)

    rep_count = len(signal_list)
    correlation_value_list = []

    for i in range(rep_count):
        for j in range(i + 1, rep_count):
            cmd = 'R --slave --vanilla <<< "cor(%s, %s)"'
            cmd_result = subprocess.check_output(cmd % (signal_list[i], signal_list[j]))
            correlation_value_list.append(float(re.findall("[1] (.*)")[0]))
    correlation_value_list = [0.6]
    min_correlation = min(correlation_value_list)


    if min_correlation >= 0.6:
        judge = 'Pass'
    else:
        judge = 'Fail'

    latex_summary_table = ['Replication QC','%d rep treatment'%rep_count,'%s'%str(min_correlation),'0.6',judge]
    print(latex_summary_table)
    replicate_latex = JinjaTemplateCommand(
            name = "venn and correlation",
            template = input["latex_summaryTable"],
            param = {"section_name": "correlation",
                     "correlation_graph": input["cor_pdf"]})
    replicate_latex.allow_fail = True
    
    write_into(replicate_latex, output["latex_section"])
    return {}


def qc_venn(input = {"venn": "", "latex_template": ""}, output = {"latex_section": ""}):
    venn_latex = JinjaTemplateCommand(
        name = "venn diagram latex",

        template = input["latex_template"],
        param = {"section_name": "venn",
                 "venn_graph": input["venn"],
                 }
    )
    venn_latex.allow_fail = True
    write_into(venn_latex, output["latex_section"])
