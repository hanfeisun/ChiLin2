import json
import re
import subprocess
from chilin2.helpers import JinjaTemplateCommand, template_dump, json_dump, json_load

def stat_cor(input={"correlation_R":"", "cor_pdf": "", "venn": "", },
             output={"json": ""}, param=None):
    # TODO: merge this into stat_venn
    """ ReplicateQC aims to describe the similarity of replicate experiment. Venn diagram and correlation plot will be used."""

    correlation_result_r_code = open(input["correlation_R"]).read()
    signal_list = re.findall(r"[pc]\d? <- (.*)$", correlation_result_r_code)

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


    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    result_dict["stat"]["judge"] = judge

    json_dump(result_dict)

def latex_cor(input, output, param):
    json_dict = json_load(input["json"])
    latex = JinjaTemplateCommand(
        name = "correlation",
        template = input["template"],
        param = {"section_name": "correlation",
                 "correlation_graph": json_dict["input"]["cor_pdf"],
                 "render_dump": output["latex"]})
    template_dump(latex)




def stat_venn(input={"venn": ""}, output={"json",""}, param=None):
    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dump(result_dict)


def latex_venn(input, output, param):
    json_dict = json_load(input["json"])
    latex = JinjaTemplateCommand(
        name = "venn diagram latex",
        template = input["template"],
        param = {"section_name": "venn",
                 "venn_graph": json_dict["input"]["venn"],
                 "render_dump": output["latex"]})
    template_dump(latex)





