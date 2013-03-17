import re
import os
import sqlite3
from pyflow.helper import fetch
def python_fastqc_parse(input, output=None, param=None):
    data = open(input).readlines()
    sequence_length = 0
    quality_dict = {}
    quality_flag = False
    for line in data:
        if re.search(r"^Sequence length", line):
            assert sequence_length == 0
            sequence_length = int(re.findall(r"^Sequence length\t(\d+)", line)[0])
        elif re.search(r"^>>Per sequence quality", line):
            assert not quality_flag
            quality_flag = True
            continue

        if re.search(r"^>>END_MODULE", line):
            quality_flag = False

        if (not line.startswith("#")) and quality_flag:
            sequence_quality = re.findall("^(\w+)\t(\w+)", line)[0]
            quality_dict[sequence_quality[0]] = float(sequence_quality[1])
    total = sum(quality_dict.values())
    n = 0
    for item in sorted(quality_dict.items(), key=lambda e: e[0], reverse=True):
        n = n + item[1]
        if n / total > 0.5:
            median = int(item[0])
            break
    return {"sequence_length": sequence_length,
            "median": median}


def python_fastqc_dist_draw(input={"db": ""},
                                 output={"rfile": "", "pdf": ""},
                                 param = {"medians": [], "ids": []},
                                 **args):

    param = fetch(param)
    db = sqlite3.connect(input["db"]).cursor()
    db.execute("select median_quality from fastqc_info")
    history_data = ",".join([str(i[0]) for i in db.fetchall()])

    with open(output["rfile"],"w") as rf:
        rf.write("""
        sequence_quality_score<-c( {locals[history_data]} )
        ecdf(sequence_quality_score)->fn
        fn(sequence_quality_score)->density
        cbind(sequence_quality_score,density)->fndd
        fndd2<-fndd[order(fndd[,1]),]
        pdf('{output[pdf]}')
        plot(fndd2,type='b',pch=18,col=2,main='Sequence Quality Score Cumulative Percentage',
              ylab='cummulative density function of all public data')
        """.format(output=output, input=input, locals=locals()))
        col = ['#FFB5C5', '#5CACEE', '#7CFC00', '#FFD700', '#8B475D', '#8E388E', '#FF6347', '#FF83FA', '#EEB422',
               '#CD7054']
        pch = [21, 22, 24, 25, 21, 22, 24, 25, 21, 22, 24, 25, 21, 22, 24, 25]
        for idx, median in enumerate(param["medians"]):
            rf.write("points({cnt},fn({cnt}),pch={pch},bg='{bg}')\n".format(
                cnt=median, pch=pch[idx], col=col[idx], bg=col[idx]
            ))

        rf.write("legend('topleft',c({ids}),pch=c({pchs}),pt.bg=c({cols})".format(
            ids=",".join(str(i) for i in param["ids"]),
            pchs=",".join(str(i) for i in pch[:len(param["ids"])]),
            cols=",".join(str(i) for i in col[:len(param["ids"])])
        ))
        rf.write("dev.off()")
    os.system('Rscript %s' % output["rfile"])


