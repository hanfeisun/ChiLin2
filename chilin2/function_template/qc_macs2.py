import os
import sqlite3
import math
from pyflow.command import ShellCommand, ThrowableShellCommand

from chilin2.jinja_template_render import JinjaTemplateCommand, write_and_run_Rscript, write_into


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
def _qc_peak_summary_parse(input={"macs2_peaks_xls": ""},
                              param={"id": ""}):
    """Basic statistic of peak calling result."""
    name = 'dataset'+ param["id"]
    with open(input["macs2_peaks_xls"],"rU") as fhd:
        fold_enrichment_list = []
        for i in fhd:
            i = i.strip()
            if i.startswith("# qvalue cutoff"):
                q_value_cutoff = float(i.split('=')[1])
                continue

            if i.startswith("#") or i.startswith("chr\t") or not i:
                continue

            fold_enrichment = float(i.split("\t")[7]) #8th column is fold change
            fold_enrichment_list.append(fold_enrichment)


        fold_greater_than_10_peaks = [x for x in fold_enrichment_list if x >= 10]
        total_peak_count = len(fold_enrichment_list)

        fold_greater_than_10_peaks_count = len(fold_greater_than_10_peaks)
    peaks_summary = [name, q_value_cutoff, total_peak_count,fold_greater_than_10_peaks_count]

#    latex_summary_table = {"desc":'Peaks number with fold change greater than 10X  ',
#                         "data": name,
#                         "value": fold_greater_than_10_peaks_count,
#                         "cutoff":1000}

    return {"peaks_summary":peaks_summary,"fold_gt_10_peaks_count":fold_greater_than_10_peaks_count}


def qc_high_confident_peaks_draw(input={"macs2_peaks_xls": "", "latex_template": "","db": "", "R_template": ""},
                             output={"rfile": "", "latex_section": "", "pdf": ""},
                             param={"id": ""}):
    """ cummulative plot of peaks fold change greater than 10"""
#    param = fetch(param)
    peaks_summary_result = _qc_peak_summary_parse(input=input, param = param)

    name = [param["id"]]
    db = sqlite3.connect(input["db"]).cursor()
    db.execute("select peak_fc_10 from peak_calling")
    historyData = [math.log(i[0]+0.001,10) for i in db.fetchall() if i[0] > 0]

    high_confident_peaks_R = JinjaTemplateCommand(name="highpeaksQC",
        template=input["R_template"],
        param={'historic_data': historyData,
               'current_data': [math.log(peaks_summary_result["fold_gt_10_peaks_count"]+0.01,10)],    #
               'ids': name,
               'cutoff': 3,
               'main': 'High confidence peaks distribution',
               'xlab': 'log(Number of Peaks fold greater than 10)',
               'ylab': 'fn(log(Number of Peaks fold greater than 10))',
               "pdf": output["pdf"]})

    write_and_run_Rscript(high_confident_peaks_R, output["rfile"])
    high_confident_latex = JinjaTemplateCommand(
        name = "high confident latex",
        template = input["latex_template"],
        param = {"section_name": "Peakcalling",
                 "peak_summary_table": peaks_summary_result["peaks_summary"],
                 "high_confident_peak_graph": output["pdf"]})
    write_into(high_confident_latex, output['latex_section'])
    return {}

def _qc_redundant_rate_parse(input={"all_peak_xls":[]}, param = {"ids": []}):
    redundant_ratio_list = []
    for a_xls in input["all_peak_xls"]:
        with open(a_xls,"rU" ) as fhd:
            for i in fhd:
                if i.startswith("# Redundant rate in treatment: "):
                    redundant_ratio_list.append(float(i.split(':')[1]))
                    break
    return redundant_ratio_list

def qc_non_redundant_rate_draw(input={"all_peak_xls": [],"db": "", "latex_template": "","R_template": ""},
                    output={"rfile": "", "latex_section": "", "pdf": ""},
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
    
    write_and_run_Rscript(redundant_rate_R, output["rfile"])
    redundant_rate_latex = JinjaTemplateCommand(
        name="redunRateQC",
        template=input["latex_template"],
        param = {"section_name": "redundant",
                 "redundant_ratio_graph": output["pdf"]})

    write_into(redundant_rate_latex, output['latex_section'])
    return {}
