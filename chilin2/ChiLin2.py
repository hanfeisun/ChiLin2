#!/usr/bin/env python3
import os
import sys
import argparse
from os.path import basename
from pyflow.command import ShellCommand, PythonCommand
from pyflow.workflow import Workflow, attach_back
from pyflow.helper import LazyCollector
from chilin2.config import ChiLinConfig
from pkg_resources import resource_filename
from chilin2.python_templates.fastqc import (
    python_fastqc_dist_draw,
    python_fastqc_parse
)
import re

class FriendlyArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit()

def parse_args(args=None):
    """
    If args is None, argparse will parse from sys.argv
    """
    description = "ChiLin :  ChIP-seq pipeline"
    parser = FriendlyArgumentParser(description = description)
    sub_parsers = parser.add_subparsers(help = "sub-command help", dest = "sub_command")

    template_parser = sub_parsers.add_parser("gen",  help = "generate a template of config file",
        description = "ChiLin-gen: A config template generator for ChiLin")
    template_parser.add_argument("-s","--species", choices = ("hg19", "mm9"), required = True)
    template_parser.add_argument("-t", dest = "atype", choices = ("Dnase", "Histone", "TF"), required = True,
        help = "the most important option for ChiLin specify the analysis type and the shiftsize {Dnase: 50, Histone and TF:73} of MACS2")

    parser_run = sub_parsers.add_parser("run", help = "run pipeline using a config file",
        description = "ChiLin-run: Run ChiLin pipeline using a config file")
    parser_run.add_argument("-c","--config", dest="config", required = True,
        help = "specify the config file to use", )
    parser_run.add_argument("-m", dest = "cor_method", default = "mean",choices = ("mean",),
        help = "specify method for correlation plot")
    parser_run.add_argument("-p", dest = "top_peaks", type = int, default = 5000,
        help = "specify peaks number for CEAS")
    parser_run.add_argument("--threads", dest = "max_threads", type = int, default = 1, choices = range(1,9),
        help = "How many threads can be used")
    parser_run.add_argument("-v", "--verbose-level", dest="verbose_level", type=int, default=2)
    parser_run.add_argument("--dry-run", dest="dry_run", action="store_true", default=False)
    parser_run.add_argument("--allow-dangling", dest="allow_dangling", action="store_true", default=False)
    parser_run.add_argument("--resume", dest="resume", action="store_true", default=False)
    parser_run.add_argument("--debug", help = "debug mode", action = "store_true", default = False)
    return parser.parse_args(args)

def step1_prepare_groom_sequencing_files(workflow, conf):

    # Return the file name pair (raw, target) that can't be groomed
    # If a file is neither a bam nor a fastq, it's returned in a list.
    # This files (maybe BED files) can be captured and processed by other tools later
    not_groomed = []
    for raw, target in conf.sample_map:
        if re.search(r"\.bam", raw, re.I):
            attach_back(workflow,
                ShellCommand(
                    "{tool} -i {input} -fq {output}",
                    tool = "bamToFastq",
                    input = raw,
                    output = target + ".fastq",
                    name = "groom"))

        elif re.search(r"\.(fastq|fq)", raw, re.I):
            attach_back(workflow,
                ShellCommand(
                    "cp {input} {output}",
                    input = raw,
                    output = target + ".fastq",
                    name = "groom"))

        else:
            print(raw, " is neither fastq nor bam file. Skip grooming.")
            not_groomed.append(raw,target)
    return not_groomed


def step2_prepare_raw_QC(workflow, conf):
    parsed_result_list = []
    for raw, target in conf.sample_map:
        fastqc_run = attach_back(workflow,
            ShellCommand(
                "{tool} {input} --extract -t {param[threads]} -o {output[target_dir]}",
                input = target + ".fastq",
                output = {"target_dir": conf.target_dir,
                          "summary": target + "_fastqc/fastqc_data.txt"},
                tool = "fastqc",
                param = {"threads": 4}))
        fastqc_run.update(param=conf.items("fastqc"))

        fastqc_parse = attach_back(workflow,
            PythonCommand(
                python_fastqc_parse,
                input = fastqc_run.output["summary"],
                param = {"id": basename(target)}))
        parsed_result_list.append(fastqc_parse)

    attach_back(workflow,
        PythonCommand(
            python_fastqc_dist_draw,
            input = {"db": resource_filename('chilin2', 'db/ChiLinQC.db')},
            output = {"rfile": conf.target_prefix + "_historical_draw.r",
                      "pdf": conf.target_prefix + "_historical_draw.pdf"},
            param = LazyCollector(
                parsed_result_list,
                medians = lambda _: _.result["median"],
                ids = lambda _: _.param["id"])))


def step3_prepare_bowtie_map(workflow, conf):
    for raw, done in conf.sample_map:
        bowtie_map = attach_back(workflow,
            ShellCommand(
                "{tool} -p {param[threads]} -S -m {param[max_align]} \
                {param[genome_index]} {input[fastq]} {output[sam]}",
                input = {"genome_dir": os.path.dirname(conf.get_path("lib", "genome_index")),
                         "fastq": done + ".fastq"},
                output = {"sam": done + ".sam"},
                tool = "bowtie",
                param = {"threads":4,
                         "max_align":1,
                         "genome_index": conf.get_path("lib", "genome_index")}))

        bowtie_map.update(param=conf.items("bowtie"))

def step4_prepare_macs2_peakcall(workflow, conf):
    for raw,done in conf.sample_map:
        attach_back(workflow,
            ShellCommand(
                "{tool} view -bt {input[chrom_len]} {input[sam]} -o {output[bam]}",
                tool="samtools",
                input={"sam":done + ".sam", "chrom_len": conf.get_path("lib", "chrom_len")},
                output={"bam":done+ ".bam"}))

    merge_bams_treat = attach_back(workflow,
        ShellCommand(
            "{tool} merge {output[merged]} {param[bams]}",
            tool = "samtools",
            input = [done+".bam" for _, done in conf.treatment_map],
            output= {"merged": conf.target_prefix + "_treatment.bam"}))

    # Here I use `merge_bams_treat` as a prototype to avoid duplication
    merge_bam_control = merge_bams_treat.clone
    merge_bam_control.input = [done+".bam" for _, done in conf.control_map]
    merge_bam_control.output = {"merged": conf.target_prefix + "_control.bam"}

    # I have to set `params` separately as it depends on `input` parameter
    merge_bams_treat.param = {"bams":" ".join(merge_bams_treat.input)}
    merge_bam_control.param = {"bams":" ".join(merge_bams_treat.input)}

    # remember to attach_back into the workflow for prototype
    attach_back(workflow, merge_bam_control)

    # Most complicated command
    macs2_on_merged = attach_back(workflow,
        ShellCommand(
            "{tool} callpeak -B -q 0.01 --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel \
            -t {input[treat]} -c {input[control]} -n {param[description]}",
            tool = "macs2",
            input = {"treat": conf.target_prefix + "_treatment.bam",
                     "control": conf.target_prefix + "_control.bam"},
            output = {"peaks": conf.target_prefix + "_peaks.bed",
                      "summit": conf.target_prefix + "_summits.bed",
                      "treat_bdg": conf.target_prefix + "_treat_pileup.bdg",
                      "ENCODE": conf.target_prefix + "_peaks.encodePeak",
                      "peaks_xls": conf.target_prefix + "_peaks.xls",
                      "control_bdg": conf.target_prefix + "_control_lambda.bdg"},
            param = {"description": conf.target_prefix,
                     "keep_dup": 1,
                     "shiftsize":73},
            name= "macs2_callpeak_merged"))
    
    macs2_on_merged.update(param=conf.items("macs2"))
    

    ## filter bdg file to remove over-border coordinates
    bdg_trim_control = attach_back(workflow,
                           ShellCommand(
            '{tool} intersect -a {input} -b {param[chrom_bed]} -wa -f 1.00 > {output}',
            tool = "bedtools",
            input = conf.target_prefix + "_control_lambda.bdg",
            output = conf.target_prefix + "_control_lambda.bdg.tmp",
            param = {'chrom_bed': conf.get("lib", "chrom_bed")},
            name = "bedGraph filtering"
            ))

    ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
    bdg_trim_treat = bdg_trim_control.clone
    bdg_trim_treat.input = conf.target_prefix + "_treat_pileup.bdg"
    bdg_trim_treat.output = conf.target_prefix + "_treat_pileup.bdg.tmp"
    attach_back(workflow, bdg_trim_treat)
    
    bdg2bw_treat = attach_back(workflow,
        ShellCommand(
            "{tool} {input[bdg]} {input[chrom_len]} {output[bw]}",
            tool="bedGraphToBigWig",
            input={"bdg": conf.target_prefix + "_control_lambda.bdg.tmp",
                   "chrom_len": conf.get("lib", "chrom_len")},
            output={"bw": conf.target_prefix + "_control.bw"},
            name= "bdg_to_bw"))

    # prototype used here to do the similar thing on treatment file
    bdg2bw_control = bdg2bw_treat.clone
    bdg2bw_control.input["bdg"] = conf.target_prefix + "_treat_pileup.bdg.tmp"
    bdg2bw_control.output["bw"] = conf.target_prefix + "_treat.bw"
    attach_back(workflow, bdg2bw_control)

def step4_prepare_macs2_peakscall_on_rep(workflow, conf):
    # Though macs command already exists, I choose not to use prototype here
    # Because the prototype definition and usage might be far from each other, making codes not readable
    for _, target in conf.treatment_map:
        macs2_on_rep = attach_back(workflow,
            ShellCommand(
                "{tool} callpeak -B -q 0.01 --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel \
                -t {input[treat]} -c {input[control]} -n {param[description]}",
                tool = "macs2",
                input = {"treat": target + ".bam",
                         "control": conf.target_prefix + "_control.bam"},
                output = {"peaks": target + "_peaks.bed",
                          "summit": target + "_summits.bed",
                          "treat_bdg": target + "_treat_pileup.bdg",
                          "ENCODE": target + "_peaks.encodePeak",
                          "peaks_xls": target + "_peaks.xls",
                          "control_bdg": target + "_control_lambda.bdg"},
                param = {"description": target, "keep_dup": 1, "shiftsize":73},
                name= "macs2_callpeak_rep"))
        macs2_on_rep.update(param=conf.items("macs2"))

        ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
        ## filter bdg file to remove over-border coordinates
        bdg_trim_controlrep = attach_back(workflow,
                                          ShellCommand(
                '{tool} intersect -a {input} -b {param[chrom_bed]} -wa -f 1.00 > {output}',
                tool = "bedtools",
                input = target + "_control_lambda.bdg",
                output = target + "_control_lambda.bdg.tmp",
                param = {'chrom_bed': conf.get("lib", "chrom_bed")},
                name = "bedGraph replicate filtering"
                ))
        bdg_trim_treatrep = bdg_trim_controlrep.clone
        bdg_trim_treatrep.input = target + "_treat_pileup.bdg"
        bdg_trim_treatrep.output = target + "_treat_pileup.bdg.tmp"
        attach_back(workflow, bdg_trim_treatrep)
        
        bdg2bw_treatrep = attach_back(workflow,
                                      ShellCommand(
                "{tool} {input} {param[chrom_len]} {output}",
                tool="bedGraphToBigWig",
                input= target + "_treat_pileup.bdg.tmp",
                output = target + "_treat.bw",
                param = {"chrom_len": conf.get("lib", "chrom_len")},name= "bdg_to_bw"))

        # prototype used here to do the similar thing on treatment file
        bdg2bw_controlrep = bdg2bw_treatrep.clone
        bdg2bw_controlrep.input = target + "_treat_pileup.bdg.tmp"
        bdg2bw_controlrep.output = target + "_treat.bw"
        attach_back(workflow, bdg2bw_controlrep)

def step4_prepare_macs2_venn_on_rep(workflow, conf):
    # awk and bedClip to remove outlier for venn and correlation plot
    for _, target in conf.treatment_map:
        bed_filter = attach_back(workflow,
                ShellCommand(
                "{tool} '{{if ($2 >= 0 && $2 < $3) print}}' {input} > {output}",
                tool = "awk",
                input = target + "_peaks.bed",
                output = target + "_peaks.bed.tmp",
                name = "filter bed files"))

        # prototype used here to do the similar thing on bedclip
        bedclip = attach_back(workflow,
                              ShellCommand(
                              template = "{tool} {input} {param[chrom_len]} {output}",
                              tool = "bedClip",
                              input  = target + "_peaks.bed.tmp",
                              output = target + "_peaks.bed",
                              param = {'chrom_len': conf.get_path("lib", "chrom_len")},
                              name = "bedclip filter"))
        bedclip.allow_fail = True


    venn_on_peaks = attach_back(workflow,
        ShellCommand(
            "{tool} -t Overlap_of_Replicates {param[beds]} && \
            mv venn_diagram.png {output}",
            tool = "venn_diagram.py",
            input = [target + "_peaks.bed.tmp" for _, target in conf.treatment_map],
            output = conf.target_prefix + "venn.png", name= "venn_diagram"))
    venn_on_peaks.param = {"beds": " ".join(venn_on_peaks.input)}
    venn_on_peaks.allow_fail = True

def step4_prepare_macs2_cor_on_rep(workflow, conf):

    # TODO (Qian): fix output.R into output[R] as output.R only works in jinja
    cor_on_bw = attach_back(workflow,
        ShellCommand(
            template =
            """{tool} \
            -s {param[wig_correlation_step]}  \
            --min-score {param[wig_correlation_min]} --max-score {param[wig_correlation_max]} \
            -r {output[R]} {param[bw]}  {param[rep]} && \
            mv {output[R]}.pdf {output[pdf]}""",
            tool = "bigwig_correlation.py",
            input = [target + "_treat.bw" for _, target in conf.treatment_map],
            output = {"R": conf.target_prefix + "_cor.R", "pdf": conf.target_prefix + "_cor.pdf"},
            param = {"wig_correlation_method": "mean",
                     "wig_correlation_min": 2,
                     "wig_correlation_max": 50},
            name = "cor_on_bw"))
    cor_on_bw.param["bw"] = " ".join(cor_on_bw.input)
    cor_on_bw.param["rep"] = " ".join([" -l replicate_%s"%(x+1) for x in range(len(conf.treatment_map))])
    cor_on_bw.update(param = conf.items("correlation"))


def step5_prepare_ceas_annotation(workflow, conf):
    get_top_peaks = attach_back(workflow,
        ShellCommand(
            "{tool} -r -g -k 5 {input} | head -n {param[peaks]} > {output}",
            tool = "sort",
            input = conf.target_prefix + "_peaks.bed",
            output = conf.target_prefix + "_peaks_top.bed",
            param = {"peaks": 3000},
            name= "top summits for conservation"))
    get_top_peaks.update(param = conf.items("ceas"))

    ceas = attach_back(workflow,
        ShellCommand(
            "{tool} -g {input[gene_table]} \
            -l {input[chrom_len]} \
            --name {param[description]} -b {input[bed]} -w {input[bw]}",
            tool="ceasBW",
            input={"bed": conf.target_prefix + "_peaks.bed",
                   "bw": conf.target_prefix + "_treat.bw",
                   "gene_table": conf.get_path("lib", "gene_table"),
                   "chrom_len": conf.get_path("lib", "chrom_len")
                   },
            output={"R": conf.target_prefix + "_ceas.R",
                    "pdf": conf.target_prefix + "_ceas.pdf"},
            param = {"description": conf.target_prefix},
            name= "ceas"))
    ceas.update(param = conf.items("ceas"))
    ceas.allow_fail = True

def step5_prepare_phast_conservation_annotation(workflow, conf):
    get_top_peaks = attach_back(workflow,
        ShellCommand(
            "{tool} -r -g -k 5 {input} | head -n {param[peaks]} > {output}",
            tool = "sort",
            input = conf.target_prefix + "_summits.bed",
            output = conf.target_prefix + "_summits_topconserv.bed",
            param = {'peaks': 5000}, name= "top summits for conservation"))
    get_top_peaks.update(param=conf.items('conservation'))

    conservation = attach_back(workflow,
        ShellCommand(
            "{tool} -t Conservation_at_summits \
            -d {input[phast]} -l Peak_summits {input[bed]} -w {param[width]}",
            tool = "conservation_plot.py",
            input = {"bed" : conf.target_prefix + "_summits_topconserv.bed",
                     "phast" : conf.get_path("lib", "phast")},
            output = ["tmp.pdf", "tmp.R"],
            param = {"width": 4000},
            name= "conservation"))
    conservation.update(param = conf.items('conservation'))

    attach_back(workflow,
        ShellCommand(
            "{tool} -resize 500x500 -density 50  {input[pdf]} {output[pdf]} && mv {input[R]} {output[R]}",
            tool = "convert",  ## width differs histone mark and TF
            input = {"pdf": "tmp.pdf", "R": "tmp.R"},
            output = {"pdf": conf.target_prefix + "conserv.pdf", "R": conf.target_prefix + "conserv.R"},
            name= "convert pdf to png",))



def step5_prepare_mdseqpos_annotation(workflow, conf):

    # This will work for only Human and Mouse
    # MDseqpos take uniform input, we need to remove random sequence part and get top 1000 summits
    attach_back(workflow,
        ShellCommand(
            '{tool} "/^chr[1-22XY]/" {input} |sort -r -g -k 5|head -n {param[peaks]} > {output}',
            tool = "awk",
            input = conf.target_prefix + "_summits.bed",
            output = conf.target_prefix + "_summits_topmdfilter.bed",
            param = {'peaks': 1000}, name= "filter summits for motif")
    ).update(param = conf.items("seqpos"))

    attach_back(workflow,
        ShellCommand(
            "{tool} -d  -w 600  -p 0.001  -m cistrome.xml  -s {param[species]} {input} {param[version]}",
            "MDSeqPos.py",
            input = conf.target_prefix + "_summits_topmdfilter.bed",
            output = "results",
            param = {"species": "hs", "version": "hg19"}, name= "motif finding"
        )).update(param = conf.items("seqpos"))
    attach_back(workflow,
        ShellCommand(
            "{tool} {input} {output}",
            "mv",
            input = "results",
            output = conf.target_prefix + "_seqpos",
            name= "mv seqpos"
        ))
def step6_prepare_report_summary(workflow, conf):
    attach_back(workflow,
        ShellCommand(name= "report"
        )
    )
        
def create_workflow(args, conf):
    """
    :type conf:ChiLinConfig
    """
    workflow = Workflow(name="Main")

    # If the output doesn't exist yet, create one
    attach_back(workflow,
        ShellCommand(
            "if [ ! -d '{output}' ]; then mkdir -p {output}; fi",
            output = conf.target_dir))

    # Whether there are replicates for treatment group
    have_reps = len(conf.treatment_map) >= 2


    step1_prepare_groom_sequencing_files(workflow, conf)
    step2_prepare_raw_QC(workflow, conf)
    step3_prepare_bowtie_map(workflow, conf)
    step4_prepare_macs2_peakcall(workflow, conf)
    if have_reps:
        step4_prepare_macs2_peakscall_on_rep(workflow, conf)
        step4_prepare_macs2_venn_on_rep(workflow, conf)
        step4_prepare_macs2_cor_on_rep(workflow, conf)
    step5_prepare_ceas_annotation(workflow, conf)
    step5_prepare_phast_conservation_annotation(workflow, conf)
    step5_prepare_mdseqpos_annotation(workflow, conf)
    return workflow

def main(args=None):
    args = parse_args(args)
    print("Arguments:", args)
    conf = ChiLinConfig(args.config)
    workflow = create_workflow(args, conf)
    workflow.set_option(
        verbose_level = args.verbose_level,
        dry_run_mode = args.dry_run,
        resume = args.resume,
        allow_dangling = args.allow_dangling)

    workflow.invoke()


if __name__ == "__main__":
    main()
