import json
import sqlite3
import math

from chilin2.helpers import JinjaTemplateCommand, template_dump, r_exec, json_load, json_dump, underline_to_space


def _peaks_parse(input):
    total = 0
    fc20n = 0
    fc10n = 0
    peaks_info = {}
    with open(input) as peaks_xls:
        for line in peaks_xls:
            if line.startswith('# tags after filtering in treatment'):
                # tags after filtering in treatment: 13438948
                peaks_info["treat_unic"] = int(line.strip().split(":")[1])
            if line.startswith("# Redundant rate in treatment: "):
                peaks_info["treat_unic_ratio"] = 1 - float(line.split(':')[1])

            if line.startswith('# tags after filtering in control'):
                peaks_info["control_unic"] = int(line.strip().split(":")[1])
            if line.startswith("# Redundant rate in control: "):
                peaks_info["control_unic_ratio"] = 1 - float(line.split(':')[1])

            if line.startswith('# d'):
                peaks_info["distance"] = int(line.strip().split("=")[1])
            if line.strip() != "" and not line.startswith("#") and not line.startswith("chr\t"):
                l = line.strip().split("\t")
                total += 1
                ## column 7th denotes fold change value
                fc = float(l[7])
                if fc >= 20:
                    fc20n += 1
                if fc >= 10:
                    fc10n += 1
            if line.startswith("# qvalue cutoff"):
                q_value_cutoff = float(line.split('=')[1])
            if line.startswith("# d"): # parse shift-size, # d =
                shift_size = int(line.strip().split("=")[1])/2

    peaks_info["totalpeak"] = total
    peaks_info["peaksge20"] = fc20n
    peaks_info["peaksge10"] = fc10n
    peaks_info["peaksge20ratio"] = peaks_info["peaksge20"] / peaks_info["totalpeak"]
    peaks_info["peaksge10ratio"] = peaks_info["peaksge10"] / peaks_info["totalpeak"]
    peaks_info["qvalue"] = q_value_cutoff
    peaks_info["shiftsize"] = shift_size
    return peaks_info

def stat_macs2(input={"macs2_peaks_xls": "", "db": "", "template": ""},
               output={"R": "", "json": "", "pdf": ""},
               param={"id": ""}):
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict["stat"] = _peaks_parse(input["macs2_peaks_xls"])
    json_dict["stat"]["cutoff"] = {"uni_loc": 5000000, "high_conf_peaks": 1000}
    json_dict["stat"]["judge"] = {"uni_loc": "Pass" if  json_dict["stat"]["treat_unic"] > 5000000 else "Fail",
                                  "high_conf_peaks": "Pass" if json_dict["stat"]["peaksge10"] >= 1000 else "Fail"}

    name = [param["id"]]
    db = sqlite3.connect(input["db"]).cursor()
    db.execute("select peak_fc_10 from peak_calling")
    historyData = [math.log(i[0] + 0.001, 10) for i in db.fetchall() if i[0] > 0]

    high_confident_peaks_r = JinjaTemplateCommand(name="highpeaksQC",
        template=input["template"],
        param={'historic_data': historyData,
               'current_data': [math.log(json_dict["stat"]["peaksge10"] + 0.01, 10)], #
               'ids': name,
               'cutoff': 3,
               'main': 'High confidence peaks distribution',
               'xlab': 'log(Number of Peaks fold greater than 10)',
               'ylab': 'fn(log(Number of Peaks fold greater than 10))',

               "pdf": output["pdf"],
               "render_dump": output["R"]})

    template_dump(high_confident_peaks_r)
    r_exec(high_confident_peaks_r)

    json_dump(json_dict)


def latex_macs2(input, output, param):
    # TODO: qian work out the peaks_summary_result part


    json_dict = json_load(input["json"])

    summary = [underline_to_space(json_dict["param"]["id"]),
               json_dict["stat"]["qvalue"],
               json_dict["stat"]["totalpeak"],
               json_dict["stat"]["peaksge10"],
               json_dict["stat"]["shiftsize"]]

    high_confident_latex = JinjaTemplateCommand(
        name = "high confident latex",
        template = input["template"],
        param = {"section_name": "high_confident_peaks",
                 "peak_summary_table": summary,
                 "high_confident_peak_graph": json_dict["output"]["pdf"],
                 "render_dump": output["latex"]})

    template_dump(high_confident_latex)



def stat_macs2_on_treats(input={"all_peak_xls": [], "db": "", "template": "", "template": ""},
                         output={"R": "", "json": "", "pdf": ""},
                         param={"ids": []}):
    """ Show redundant  ratio of the dataset in all historic data
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}

    treat_unic_loc_rates = []
    treat_unic_locs = []
    control_unic_loc_rates = []
    for i in input["all_peak_xls"]:
        parsed = _peaks_parse(i)
        treat_unic_loc_rates.append ( parsed["treat_unic_ratio"])
        treat_unic_locs.append(parsed["treat_unic"])



    for id, unic_loc_rate, unic_loc in zip(param["ids"], treat_unic_loc_rates, treat_unic_locs):
        json_dict["stat"][id] = {"unic_loc_rate": unic_loc_rate,"unic_loc": unic_loc }
        json_dict["stat"][id]["cutoff"] = {"unic_loc_rate": 0.8, "unic_loc": 5000000} # non redundant rate cutoff
        json_dict["stat"][id]["judge"] = {"unic_loc_rate": "Pass" if  unic_loc_rate >= 0.8 else "Fail",
                                          "unic_loc": "Pass" if unic_loc >= 5000000 else "Fail"}

    db = sqlite3.connect(input["db"]).cursor()

    # TODO: column name redundant_rate => non-redundant rate
    db.execute("select redundant_rate from peak_calling")
    redundant_history = db.fetchall()
    historyData = [1 - float(i[0]) for i in redundant_history if i[0] != "null"]

    redundant_rate_R = JinjaTemplateCommand(name="redunRateQC",
        template=input["template"],
        param={'historic_data': historyData,
               'current_data': treat_unic_loc_rates,
               'ids': param["ids"],
               'cutoff': 0.8,
               'main': 'Non-Redundant rate',
               'xlab': 'Non-Redundant rate',
               'ylab': 'fn(Non-Redundant rate)',
               "pdf": output["pdf"],
               "render_dump": output["R"],
               "need_smooth_curve": True})

    template_dump(redundant_rate_R)
    r_exec(redundant_rate_R)

    with open(output["json"], "w") as f:
        json.dump(json_dict, f, indent=4)

def latex_macs2_on_sample(input, output, param):
    json_dict = json_load(input["json"])
    latex = JinjaTemplateCommand(
        name="redunRateQC",
        template=input["template"],
        param = {"section_name": "redundant",
                 "redundant_ratio_graph": json_dict["output"]["pdf"],
                 "render_dump": output["latex"]})

    template_dump(latex)



def stat_velcro(input={"macs2_peaks_xls": "", "velcro_peaks": ""}, output={"json": ""},
                param={}):
    peaks_info = _peaks_parse(input["macs2_peaks_xls"])
    peaks_info["velcro"] = len(open(input["velcro_peaks"], 'r').readlines())
    peaks_info['nonvelcropercentage'] = 1 - float(peaks_info["velcro"]) / peaks_info["totalpeak"]

    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    result_dict["stat"] = peaks_info
    result_dict["stat"]["cutoff"] = 0.9
    result_dict['stat']["judge"] = "Pass" if (peaks_info['nonvelcropercentage']) >= 0.9 else "Fail"
    with open(output["json"], "w") as f:
        json.dump(result_dict, f, indent=4)


def stat_dhs(input={"macs2_peaks_xls": "", "dhs_peaks": ""}, output={"json": ""},
             param={}):
    peaks_info = _peaks_parse(input["macs2_peaks_xls"])
    peaks_info["dhs"] = len(open(input["dhs_peaks"], 'r').readlines())
    peaks_info['dhspercentage'] = peaks_info["dhs"] / peaks_info["totalpeak"]

    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    result_dict["stat"] = peaks_info
    result_dict["stat"]["cutoff"] = 0.8
    result_dict["stat"]["judge"] = "Pass" if result_dict["stat"]["dhspercentage"] >= 0.8 else "Fail"
    with open(output["json"], "w") as f:
        json.dump(result_dict, f, indent=4)


