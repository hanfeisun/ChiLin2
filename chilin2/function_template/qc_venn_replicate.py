import re
import os
import subprocess

def replicate_info(input={"correlation_R":""},
                    output={"latex_SummaryTable":""}):
    """ ReplicateQC aims to describe the similarity of replicate experiment. Venn diagram and correlation plot will be used."""


    rscript_calculate_correlation = ""

    correlation_result_R_code = open(input["correlation_R"]).read()
    signal_list = re.findall(r"[pc]\d? <- (.*)$", correlation_result_R_code)

    rep_count = len(signal_list)

    correlation_value_list = []

    for i in range(rep_count):
        for j in range(i + 1, rep_count):
            cmd = 'R --slave --vanilla <<< "cor(%s, %s)"'
            cmd_result = subprocess.check_output(cmd % (signal_list[i], signal_list[j]))
            correlation_value_list.append(float(re.findall("[1] (.*)")[0]))

    min_correlation = min(correlation_value_list)


    if min_correlation >= 0.6:
        judge = 'Pass'
    else:
        judge = 'Fail'

    latex_summary_table = ['Replication QC','%s rep treatment'%rep_index,'%s'%str(cor),'0.6',judge]
    with open(output["latex_SummaryTable"],"a+") as f:
        f.write('\t'.join(latex_summary_table)+'\n')
    return {}
