import re
import os

def redraw_ceas_graph(input={"peaksxls": "","ceas_rscript":""},
                output={"rfile": "", "peakheight_and_pie_pdf": "","metagene_dist_pdf":""}):

    
    """ Describe peaks' distribution and relative position.
    input: the xls file and
    output: rfile and pdf
    """
    
    list_fold_change = []
    with open(input["peakxls"]) as xls_file:
        for line in xls_file:
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            if line.startswith("chr\t"):
                line_splited = line.split("\t")
                fold_change = float(line_splited[7]) # The 8th column is fold change
                list_fold_change.append(fold_change)

    foldchange_string = ','.join(map(str, list_fold_change))

    ceas_rscript_content = open(input["ceas_rscript"]).read()
        
    # TODO : change the regex here
    metagene_dist_r_codes = re.findall(
        r'''layout\(matrix\(c\(1, 2, 3, 3, 4, 5\), 3, 2, byrow = TRUE\),widths=c\(1, 1\),heights=c\(1, 1, 1\)\)''' +
        r'''(.*abline\(v=3000\.000000,lty=2,col=c\("black"\)\))''',
        ceas_rscript_content,
        re.DOTALL)[0] # return three figures for metagene distribution

    pie_chart_r_codes = re.findall(
        "(" +
        re.escape("par(mfcol=c(2, 2),mar=c(3, 3, 4, 2.7999999999999998),oma=c(4, 2, 4, 2))") +
        ".*)" +
        re.escape("# ChIP regions over the genome"),      # stop here

        ceas_rscript_content,
        re.DOTALL,)[0]


    with open(output["rfile"],'w') as r_file:
        
        # R codes for Peak height distribution & Pie chart
        r_file.write("pdf('%s',height=11.5,width=8.5)\n" %output["peakheight_and_pie_pdf"] )
        r_file.write('nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow=TRUE),respect=TRUE)\n')
        r_file.write('peaks_fc <- c(%s)\n' %foldchange_string)
        r_file.write('fn <- ecdf(peaks_fc)\n')
        r_file.write('density <- fn(peaks_fc)\n')
        r_file.write('fdd <- cbind(peaks_fc,density)\n')
        r_file.write('fdd1 <- fdd[order(fdd[,1],decreasing = TRUE),]\n')
        r_file.write('fdd2 <- cbind(fdd1[,1],1-fdd1[,2])\n')
        r_file.write('ma <- max(fdd1[,1])\n')
        r_file.write('mi <- min(fdd1[,1])\n')
        r_file.write("plot(fdd2,type='p',col='blue',pch=18,main='Peaks distribution',xlab='Fold change of peaks',ylab='Fn(fold change of peaks)')\n")
        r_file.write(pie_chart_r_codes)
        r_file.write('dev.off()\n')
        
        # R codes for MetaGene Distribution
        r_file.write('# the graph of MetaGene Distribution \n\n')
        r_file.write("pdf('%s',height=7,width=5.5)\n" %output["metagene_dist_pdf"])
        r_file.write("nf <- layout(matrix(c(1,2,3,3), 2, 2, byrow=TRUE), width= c(1,1),height=c(1,1),respect=TRUE)\n")
        r_file.write(metagene_dist_r_codes)
        r_file.write('dev.off()\n')

    os.system('Rscript %s' % output["rfile"])
    return {'AnnotationQC_check' : True,
            'ceas_check' : True,
            "peakheight_and_pie_pdf":output["peakheight_and_pie_pdf"],
            "metagene_dist_pdf":output["metagene_dist_pdf"]}

