import os
import sqlite3
from pyflow.command import ShellCommand, ThrowableShellCommand

from chilin2.jinja_template_render import JinjaTemplateCommand


# TODO: Check whether _check can be separated
#def _check(input,output=[]):
#    """ Check whether the quality of the dataset is ok.
#    vcb is the callback function for value"""
#    ord = lambda x:[x["desc"], x["data"], x["value"], x["cutoff"], x["test"]]
#    if len(input)!=0:
#        for i in input:
#            value = float(i["value"]),3
#            cutoff = float(i["cutoff"]),3
#            if value >= cutoff:
#                i["test"] = 'Pass'
#                output.append(ord(i))
#            else:
#                i["test"] = 'Fail'
#                output.append(ord(i))
#        return output
#
#
def _peak_summary_info(input={"macs2_xls_for_each_sample":"","name":""}, output={"latex_SummaryTable":""}, param={"id":""}):
    """Basic statistic of peak calling result."""
    name = 'dataset'+ param["id"]
    with open(input["macs2_xls_for_each_sample"],"rU" ) as fhd:
        fold_enrichment_list = []
        q_value_cutoff = "unknown"
        for i in fhd:
            i = i.strip()
            if i.startswith("# qvalue cutoff"):
                q_value_cutoff = float(i.split('=')[1])

            if i.startswith("#") or i.startswith("chr\t"):
                continue

            fold_enrichment = float(i.split("\t")[7]) #8th column is fold change
            fold_enrichment_list.append(fold_enrichment)


        fold_greater_than_10_peaks = [x for x in fold_enrichment_list if x >= 10]
        total_peak_count = len(fold_enrichment_list)

        fold_greater_than_10_peaks_count = len(fold_greater_than_10_peaks)
    peaks_summary = [name, q_value_cutoff, total_peak_count,fold_greater_than_10_peaks_count]

    latex_summary_table = {"desc":'Peaks number with fold change greater than 10X  ',
                         "data": name,
                         "value": fold_greater_than_10_peaks_count,
                         "cutoff":1000}

    latex_summary_table = _check(latex_summary_table,output)
    with open(output["latex_SummaryTable"],"a+") as f:
        f.write('\t'.join(latex_summary_table)+'\n')
    return {"peak_summary_table":peaks_summary,"fold10":fold_greater_than_10_peaks_count }

#def high_confidentPeaks_info(input={"macs2_xls_for_each_sample":"","db": ""},
#                            output={"rfile": "", "pdf": "","latex_SummaryTable":""},
#                            param = {"ids": []},
#                            **args):
#    """
#    cummulative percentage of peaks foldchange great than 10
#    """
#    param = fetch(param)
#    peaks_summary_result = _peak_summary_info(input,output)
#    name = param["ids"]
#    db = sqlite3.connect(input["db"]).cursor()
#    db.execute("select peak_fc_10 from peak_calling_tb")
#    highpeaks_history = db.fetchall()
#    historyData = [str(math.log(i[0]+0.001,10)) for i in highpeaks_history if i[0] > 0]
#    historyData = ','.join(historyData)
#
#    lg_10 = round(math.log(peaks_summary_result["fold10"],10),3)
#    Rrender = {'histroyData':historyData,
#                'value':lg_10,
#                'name':"'ratio of fold change upper than 10'",
#                'cutoff':3,
#                'pdfName':output["pdf"],
#                'main':'High confidence peaks distribution',
#                'xlab':'log(Number of Peaks fold upper than 10)',
#                'ylab':'fn(log(Number of Peaks fold upper than 10))',
#                'other':True}
#    highpeaksR = JinjaCommand(name="highpeaksQC", template="Rtemplate.tex", Rrender)
#    if highpeaksR.invoke():
#        content = highpeaksR.result().replace('%','\\%')
#    with open(output["rfile"], 'w') as f:
#        f.write(content)
#    os.system('Rscript %s' % output["rfile"])
#    latex_summary_table = {"desc":'Overlap with DHSs  ',
#                         "data":'%s'%name,
#                         "value":'%f'%dhs_ratio,
#                         "cutoff":0.8}
#
#    latex_summary_table = _check(latex_summary_table,output)
#    with open(output["latex_SummaryTable"],"a+") as f:
#        f.write('\t'.join(latex_summary_table)+'\n')
#    return {'PeakcallingQC_check':True,'peak_summary_table':peaks_summary_result["peaks_summary"],
#            'high_confident_peak_graph':output["pdf"]}
#
#



def _qc_redundant_rate_parse(input={"all_peak_xls":[]}, param = {"ids": []}):
    redundant_ratio_list = []
    for a_xls in input["all_peak_xls"]:
        with open(a_xls,"rU" ) as fhd:
            for i in fhd:
                if i.startswith("# Redundant rate in treatment: "):
                    redundant_ratio_list.append(float(i.split(':')[1]))
                    break
    return redundant_ratio_list


def qc_non_redundant_rate_draw(input={"all_peak_xls": [],"db": "", "R_template": ""},
                    output={"rfile": "", "pdf": ""},
                    param = {"ids": []}):
    """ Show redundant  ratio of the dataset in all historic data
    """

    non_redundant_rate_list = [1 - i for i in _qc_redundant_rate_parse(input=input, param = param)]

    db = sqlite3.connect(input["db"]).cursor()

    # TODO: column name redundant_rate => non-redundant rate
    db.execute("select redundant_rate from peak_calling")
    redundant_history = db.fetchall()
    historyData = [1-float(i[0]) for i in redundant_history if i[0] != "null"]

    redundant_rate_R = JinjaTemplateCommand(name="redunRateQC",
        template=input["R_template"],
        param={'historic_data': historyData,
               'current_data': non_redundant_rate_list,
               'ids': param["ids"],
               'cutoff': 0.8,
               'main': 'Non-Redundant rate',
               'xlab': 'Non-Redundant rate',
               'ylab': 'fn(Non-Redundant rate)',
               "pdf": output["pdf"],
               "need_smooth_curve": True})

    redundant_rate_R.invoke()
    with open(output["rfile"], 'w') as f:
        f.write(redundant_rate_R.result)

    ThrowableShellCommand(template = 'Rscript {input}', input= output["rfile"]).invoke()
    return {}