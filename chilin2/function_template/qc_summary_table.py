import os
from chilin2.helpers import json_load, JinjaTemplateCommand, template_dump, underline_to_space

__author__ = 'ad9075'
def _items(json_path):
    return json_load(json_path)["stat"].items()

def _s(json_path):
    return json_load(json_path)["stat"]

def _2space(id):
    return underline_to_space(id)
def _cons_summary_table(conf):
    table = []
    has_run = os.path.exists
    pre = conf.json_prefix
    js = pre + "_fastqc.json"
    if has_run(js):
        for id, stat in _items(js):
            table.append(["FastQC",
                          _2space(id), stat["median"], stat["cutoff"], stat["judge"]])

    js = pre + "_bowtie.json"
    if has_run(js):
        for id, stat in _items(js):
            table.append(["Unique mappable reads",
                          _2space(id), stat["mappable_reads"], stat["cutoff"], stat["judge"]])

    js = pre + "_macs2.json"
    if has_run(js):
        stat = _s(js)
        table.append(
            ["High confident peaks", _2space(conf.id), stat["peaksge10"], stat["cutoff"]["high_conf_peaks"],
             stat["judge"]["high_conf_peaks"]])

    js = pre + "_macs2_on_sample.json"
    if has_run(js):
        macs2s_stat = _s(js)
        for id, stat in macs2s_stat.items():
            table.append(["Unique location rate",
                          _2space(id), stat["unic_loc_rate"], stat["cutoff"]["unic_loc_rate"],
                          stat["judge"]["unic_loc_rate"]])
            table.append(["Unique location",
                          _2space(id), stat["unic_loc"], stat["cutoff"]["unic_loc"],
                          stat["judge"]["unic_loc"]])

    js = pre + "_dhs.json"
    if has_run(js):
        stat = _s(js)
        table.append(["DHS ratio",
                      _2space(conf.id), round(stat["dhspercentage"], 3), stat["cutoff"], stat["judge"]])

#    js = pre + "_cor.json"
#    if has_run(js):
#        stat = _s(js)
#        table.append(["Replicates Correlation",
#                      conf.id, stat["min_cor"], stat["cutoff"], stat["judge"]])

    js = pre + "_velcro.json"
    if has_run(js):
        stat = _s(js)
        table.append(
            ["non Velro ratio",
             _2space(conf.id), round(stat["nonvelcropercentage"], 3), stat["cutoff"], stat["judge"]])

    js = pre + "_conserv.json"
    if has_run(js):
        stat = _s(js)
        table.append(["Conservation QC",
                      _2space(conf.id), "", "", stat["judge"]])
    return table

def latex_summary_table(input, output, param):
    summary_table = _cons_summary_table(param["conf"])

    latex = JinjaTemplateCommand(
        template=input["template"],
        param={"SummaryQC" : True,
               "summary_table": summary_table,
               "render_dump": output["latex"]})
    template_dump(latex)