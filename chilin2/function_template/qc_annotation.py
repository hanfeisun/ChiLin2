import sqlite3
import math
from pyflow.command import ShellCommand, ThrowableShellCommand



def qc_velcro_ratio_draw(input={"macs2_peaks_bed": "", "velcro_peaks_bed":"","db": ""},
                      output={"rfile": "", "pdf": ""},
                      param={"id": "", "latex_SummaryTable": ""}):
    """verlcro ratio is used to describe whether the peak is credible , The lower the result is more convenience.
     The cumulative percentage plot can reflect the particularly dataset's verlcro ratio quality of all historic data."""

    db = sqlite3.connect(input["db"]).cursor()
    db.execute("select velcro_rate from peak_calling")
    historyData = [str(i[0]) for i in db.fetchall()]
    historyData = [1-float(i) for i in historyData  if i!='null']

    pointText = str(velcro_ratio*100)+'%'

    high_confident_peaks_R = JinjaTemplateCommand(name="velcroRatioQC",
        template=input["R_template"],
        param={'historic_data': historyData,
               'current_data': non_velcro_rate,    #
               'ids': ['ratio overlap with non-verlcro:%s' %pointText ],
               'cutoff': 0.9,
               'main': 'Non-velcro ratio',
               'xlab': 'Non-velcro ratio',
               'ylab': 'fn(Non-velcro ratio)',
               "pdf": output["pdf"]})

    high_confident_peaks_R.invoke()
    with open(output["rfile"], 'w') as f:
        f.write(high_confident_peaks_R.result)

    ThrowableShellCommand(template = 'Rscript {input}', input= output["rfile"]).invoke()

    return {}


def qc_DHS_rate_draw(input={"macs2_peaks_bed": "", "DHS_peaks_bed":"","db": ""},
                     output={"rfile": "", "pdf": "", "latex_SummaryTable": ""},
                     param={"id": ""}):
    """ DHS ratio indicate the percentage of peaks overlap with DHSs site.
    The function can describe  the particularly dataset's DHS ratio quality of all historic data.
    """

    name = param["id"]
    db = sqlite3.connect(input["db"]).cursor()
    db.execute("select DHS_rate from peak_calling ")
    historyData = [str(i[0]) for i in db.fetchall()]
    historyData = [i for i in historyData if i!='null']



    pointText = str(DHS_rate*100)+'%'

    high_confident_peaks_R = JinjaTemplateCommand(name="DHSRateQC",
        template=input["R_template"],
        param={'historic_data': historyData,
               'current_data': DHS_rate,
               'ids': ["'ratio overlap with DHSs:%s'"%pointText ,],
               'cutoff': 0.8,
               'main': 'Overlapped_with_DHSs',
               'xlab': 'Overlapped_with_DHSs',
               'ylab': 'fn(Overlapped_with_DHSs)',
               "pdf": output["pdf"]})

    high_confident_peaks_R.invoke()
    with open(output["rfile"], 'w') as f:
        f.write(high_confident_peaks_R.result)

    ThrowableShellCommand(template = 'Rscript {input}', input= output["rfile"]).invoke()

    latex_summary_table = {"desc":'Overlap with DHSs  ',
                         "data":'%s'%name,
                         "value":'%f'%dhs_ratio,
                         "cutoff":0.8}

    return {}

