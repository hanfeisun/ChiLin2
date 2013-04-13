#def velcro_ratio_info(input={"db": ""},
#                        output={"R": "", "pdf": ""},
#                        param = {"velcro_ratio": [], "ids": [],"latex_SummaryTable":""},
#                        **args):
#    """verlcro ratio is used to describe whether the peak is credible , The lower the result is more convenience.
#     The cumulative percentage plot can reflect the particularly dataset's verlcro ratio quality of all historic data."""
#    param = fetch(param)
#    name = param["ids"]
#    db = sqlite3.connect(input["db"]).cursor()
#    db.execute("select velcro_rate from peak_calling_tb")
#    velcro_history = db.fetchall()
#    historyData = [str(i[0]) for i in velcro_history]
#    historyData = [str(1-float(i)) for i in historyData if i!='null']
#    historyData = ','.join(historyData)
#
#    pointText = str(velcro_ratio*100)+'%'
#    Rrender = {'histroyData':historyData,
#                'value':param["velcro_ratio"],
#                'name':"'ratio overlap with non-verlcro:%s'"%pointText ,
#                'cutoff':0.8,
#                'pdfName':output["pdf"],
#                'main':'non-velcro ratio',
#                'xlab':'non-velcro ratio',
#                'ylab':'fn(non-velcro ratio)',
#                'other':True}
#    velcroRatioR = JinjaCommand(name="velcroRatioQC", template="Rtemplate.tex", Rrender)
#    if velcroRatioR.invoke():
#        content = velcroRatioR.result().replace('%','\\%')
#    with open(output["R"], 'w') as f:
#        f.write(content)
#    os.system('Rscript %s' % output["R"])
#    latex_summary_table = {"desc":'Overlap with non-velcro ',
#                         "data":name,
#                         "value":velcro_ratio,
#                         "cutoff":0.9}
#    latex_summary_table = _check(latex_summary_table,output)
#    with open(output["latex_SummaryTable"],"a+") as f:
#        f.write('\t'.join(latex_summary_table)+'\n')
#    return {'verlcro_check':True,'velcro_ratio_graph':output['pdf']}
#
#def DHS_ratio_info(self,input={"db": ""},
#                        output={"R": "", "pdf": "","latex_SummaryTable":""},
#                        param = {"dhs_ratio": [], "ids": []},
#                        **args):
#    """ DHS ratio indicate the percentage of peaks overlap with DHSs site.
#    The function can describe  the particularly dataset's DHS ratio quality of all historic data.
#    """
#    param = fetch(param)
#    name = ids[1]
#    db = sqlite3.connect(input["db"]).cursor()
#    db.execute("select DHS_rate from peak_calling_tb ")
#    dhs_history = self.db.fetchall()
#    historyData = [str(i[0]) for i in dhs_history]
#    historyData = [i for i in historyData if i!='null']
#    historyData = ','.join(historyData)
#
#    pointText = str(dhs_ratio*100)+'%'
#    Rrender = {'histroyData':historyData,
#                'value':dhs_ratio,
#                'name':"'ratio overlap with DHSs:%s'"%pointText ,
#                'cutoff':0.8,
#                'pdfName':output["pdf"],
#                'main':'Overlapped_with_DHSs',
#                'xlab':'Overlapped_with_DHSs',
#                'ylab':'fn(Overlapped_with_DHSs)',
#                'other':True}
#    DHSRatioR = JinjaCommand(name="DHSRatioQC", template="Rtemplate.tex", Rrender)
#    if DHSRatioR.invoke():
#        content = DHSRatioR.result().replace('%','\\%')
#    with open(output["rfilr"], 'w') as f:
#        f.write(content)
#
#    os.system('Rscript %s' % output["R"])
#    latex_summary_table = {"desc":'Overlap with DHSs  ',
#                         "data":'%s'%name,
#                         "value":'%f'%dhs_ratio,
#                         "cutoff":0.8}
#    latex_summary_table = _check(latex_summary_table,output)
#    with open(output["latex_SummaryTable"],"a+") as f:
#        f.write('\t'.join(latex_summary_table)+'\n')
#    return {'DHS_ratio_graph':output["pdf"]}
