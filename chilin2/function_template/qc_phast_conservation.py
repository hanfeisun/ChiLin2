import re
import os
import math
from samflow.helper import fetch
from chilin2.jinja_template_render import JinjaTemplateCommand, write_into

def _euclidian_distance(x,y):
    assert len(x) == len(y)
    sum_of_square = 0
    for x_i, y_i in zip(x,y):
        sum_of_square += pow((x_i - y_i), 2)
    distance=round(math.sqrt(sum_of_square), 4)
    return distance

# TODO (Shenglin): Why should conservationPDF use here as input?
def qc_conservation_draw(input={"conservationR": "","historical_conservation_cluster_text": "", "latex_template":""},
                        output={"rfile": "","compare_pdf": "", "pdf": "", "latex_SummaryTable":""},
                        param = {"id":""}):
    """
    For TFcenters data 1,2,3 pass, 4,5,6 fail
    For Histone center data 1,2,3,4 pass, 5,6,7,8 fail.
    """
    with open(input["conservationR"]) as conservation_r_file:
        for line in conservation_r_file:
            if re.findall(r'y0<-\S*\)',line):

                # TODO: Use more friendly regex
                value = re.findall(r'y0<-\S*\)',line)[0][6:-1]
                value = value.split(',')
            elif re.findall(r'x<-c\S*\)',line):
                xlab = line
        value = [float(i) for i in value]
        sumvalue = sum(value)
        value = [i/sumvalue for i in value]

    with open(input["historical_conservation_cluster_text"]) as historic_data:
            historyData = historic_data.readlines()

    less_than_this_id_means_pass_qc = len(historyData)/2


    scoreList = []
    for i in range(len(historyData)):
        temp = historyData[i].strip()
        line = temp.split(' ')
        line = [float(j) for j in line]
        score = _euclidian_distance(value,line)
        scoreList.append(score)
    mindist = scoreList.index(min(scoreList)) # return the index of minimum distance group

    judgevalue = historyData[mindist].strip().split(' ')
    judgevalue = [str(i) for i in judgevalue]
    value = [str(i) for i in value]
    ymax = max(value+judgevalue)
    ymin = min(value+judgevalue)
    with open(output["rfile"],'w') as f:
        f.write("pdf('%s',height=8.5,width=8.5)\n" % output["compare_pdf"])
        f.write("%s\n" % xlab)
        f.write("y1<-c(%s)\n" % ','.join(judgevalue))
        f.write("y2<-c(%s)\n" % ','.join(value))
        f.write("ymax<-(%s)\n" %ymax)
        f.write("ymin<-(%s)\n" %ymin)
        f.write("yquart <- (ymax-ymin)/4\n")
        f.write("plot(x,y2,type='l',col='red',main='Normalized Conservation_at_summits',xlab='Distance from the Center (bp)',ylab='Normalized Average Phastcons',ylim=c(ymin-yquart,ymax+yquart))\n")
        f.write("points(x,y1,type='l',col='blue',lty=2)\n")
        f.write("legend('topleft',c('original group','compared group'),lty=c(1,2),col=c('red', 'blue'))\n")
        f.write("dev.off()\n")
    os.system('Rscript %s' % output["rfile"])
    if mindist <=less_than_this_id_means_pass_qc:
        judge = 'pass'
    else:
        judge = 'fail'
    # TODO (Shenglin): conf.basis.id => input["id"]
    latex_summary_table = ['Conservation QC','dataset%s'%param["id"], '%f' % round(min(scoreList),3),' ','%s'%judge]

    # TODO (Shenglin): don't use append, use a separate file for each qc function in the same folder
    with open(output["latex_section"],"a+") as f:
        f.write('\t'.join(latex_summary_table)+'\n')

    conservation_latex = JinjaTemplateCommand(
        name = "conservation",
        template = input["latex_template"],
        param={"section_name": "conservation",
               "conservation_compare_graph": output["compare_pdf"],
               "conservation_graph": output["pdf"]})
    
    write_into(conservation_latex, output["latex_section"])
    return {}
