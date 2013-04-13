#!/usr/bin/env python3
from glob import glob

import os
import random
import subprocess
import sys
import re
import argparse

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from pkg_resources import resource_filename

from chilin2.config import ChiLinConfig
from chilin2.function_template.qc_bowtie import stat_bowtie, latex_bowtie
from chilin2.function_template.qc_ceas import stat_ceas, latex_ceas
from chilin2.function_template.qc_contamination import write_random_records, library_summary
from chilin2.function_template.qc_fastqc import stat_fastqc, latex_fastqc
from chilin2.function_template.qc_macs2 import stat_macs2, stat_velcro, stat_dhs, stat_macs2_on_sample, latex_macs2, latex_macs2_on_sample
from chilin2.function_template.qc_mdseqpos import stat_seqpos, latex_seqpos
from chilin2.function_template.qc_conservation import stat_conservation, latex_conservation
from chilin2.function_template.qc_venn import stat_venn, stat_cor, latex_venn, latex_cor
from chilin2.helpers import latex_end


ChiLinQC_db = resource_filename("chilin2", "db/ChiLinQC.db")
r_template = resource_filename("chilin2", "jinja_template/R_culmulative_plot.R.jinja2")
latex_template = resource_filename("chilin2", "jinja_template/Latex_summary_report.jinja2")
plain_template = resource_filename("chilin2", "jinja_template/text_summary_report.jinja2")


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
    parser = FriendlyArgumentParser(description=description)
    sub_parsers = parser.add_subparsers(help="sub-command help", dest="sub_command")

    template_parser = sub_parsers.add_parser("gen", help="generate a template of config file",
        description="ChiLin-gen: A config template generator for ChiLin")
    template_parser.add_argument("-s", "--species", choices=("hg19", "mm9"), required=True)
    template_parser.add_argument("-t", dest="atype", choices=("Dnase", "Histone", "TF"), required=True,
        help="the most important option for ChiLin specify the analysis type and the shiftsize {Dnase: 50, Histone and TF:73} of MACS2")

    parser_run = sub_parsers.add_parser("run", help="run pipeline using a config file",
        description="ChiLin-run: Run ChiLin pipeline using a config file")
    parser_run.add_argument("-c", "--config", dest="config", required=True,
        help="specify the config file to use", )

    parser_run.add_argument("--from", dest="start_step", default=0, type=int,
        help="Only step after this number will be processed")
    parser_run.add_argument("--to", dest="end_step", default=100, type=int,
        help="Only step before this number will be processed ")
    parser_run.add_argument("--skip", dest="skip_step", default="",
        help="Steps to skip, use comma as seperator")

    parser_run.add_argument("-v", "--verbose-level", dest="verbose_level", type=int, default=2)
    parser_run.add_argument("--dry-run", dest="dry_run", action="store_true", default=False)
    parser_run.add_argument("--allow-dangling", dest="allow_dangling", action="store_true", default=False)
    parser_run.add_argument("--resume", dest="resume", action="store_true", default=False)
    parser_run.add_argument("--debug", help="debug mode", action="store_true", default=False)
    parser_run.add_argument("--remove", dest="clean", action="store_true", default=False)
    return parser.parse_args(args)


def make_copy_command(orig, dest):
    return ShellCommand(
        "cp {input} {output}",
        input=orig,
        output=dest,
        name="copy")


def _groom_sequencing_files(workflow, conf):
    """
    Return the file name pair (raw, target) that can't be groomed (neither fastq nor bam)
    These files (maybe BED files) can be captured and processed by other tools later
    """
    not_groomed = []
    for raw, target in conf.sample_pairs:
        if re.search(r"\.bam", raw, re.I):
            attach_back(workflow,
                ShellCommand(
                    "{tool} -i {input} -fq {output}",
                    tool="bamToFastq",
                    input=raw,
                    output=target + ".fastq",
                    name="groom"))

        elif re.search(r"\.(fastq|fq)", raw, re.I):
            attach_back(workflow, make_copy_command(orig=raw, dest=target + ".fastq"))
        else:
            print(raw, " is neither fastq nor bam file. Skip grooming.")
            not_groomed.append([raw, target])

def prepare_library_contamination(workflow, conf):
    """
    :param workflow: main stream workflow
    :param conf: conf parsed from config
    :return: none
    estimate library contamination ratios
    """
    ## bowtie mapping back to
    for target in conf.sample_targets:
        random_reads = attach_back(workflow,
            PythonCommand(write_random_records,
            input = {"fastq": target + ".fastq"},
            output = {"fastq_sample": target + ".fastq.subset"},
            param = {"random_number": "100000"}))     ## use 100kb random reads

        for species in dict(conf.items("contaminate_genomes")):
            attach_back(workflow,
                ShellCommand(
                "{tool} -p {param[threads]} -S -m {param[max_align]} \
                {param[genome_index]} {input[fastq]} {output[sam]} 2> {output[bowtie_summary]}",
                input={"genome_dir": os.path.dirname(conf.get_path("lib", "genome_index")),
                       "fastq": target + ".fastq.subset"},
                output={"sam": target + ".sam.subset",
                        "bowtie_summary": target + species + "_subset_bowtie_summary.txt"},
                tool="bowtie",
                param={"threads": 4,
                       "max_align": 1,
                       "genome_index": conf.get_path("contaminate_genomes", species)}))

def prepare_library_contaminationTex(workflow, conf):
    all_contaminate_files = []
    for target in conf.sample_targets:
        all_contaminate_files.append([ target + species + "_subset_bowtie_summary.txt" for species in dict(conf.items("contaminate_genomes")) ])
    Tex_step = attach_back(workflow,
        PythonCommand(library_summary,
             input = {"latex_template": latex_template},
             output = {"latex_section": conf.prefix + "_library_contamination_summary"},
             param = {"samples": conf.sample_bases,
                      "bowtie_summary": all_contaminate_files,
                      "id": conf.id,
                      "species": [os.path.basename(conf.get("contaminate_genomes", species)) for species in dict(conf.items("contaminate_genomes"))] }))
    return Tex_step.output["latex_section"]

def _raw_QC(workflow, conf):
    for target in conf.sample_targets:
        fastqc_run = attach_back(workflow,
            ShellCommand(
                "{tool} {input} --extract -t {param[threads]} -o {output[target_dir]}",
                input=target + ".fastq",
                output={"target_dir": conf.target_dir,
                        "fastqc_summary": target + "_fastqc/fastqc_data.txt"},
                tool="fastqc",
                param={"threads": 4}))
        fastqc_run.update(param=conf.items("fastqc"))

    attach_back(workflow,
        PythonCommand(
            stat_fastqc,
            input={"db": ChiLinQC_db,
                   "fastqc_summaries": [target + "_fastqc/fastqc_data.txt" for target in conf.sample_targets],
                   "template": r_template},
            output={"R": conf.prefix + "_raw_sequence_qc.R",
                    "pdf": conf.prefix + "_raw_sequence_qc.pdf",
                    "json": conf.json_prefix + "_fastqc.json"},
            param={"ids": conf.sample_bases,
                   "id": conf.id}))


def _raw_QC_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_fastqc,
            input={"json": conf.json_prefix + "_fastqc.json",
                   "template": latex_template},
            output={"latex": conf.latex_prefix + "_fastqc.latex"}))


def _bowtie(workflow, conf):
    for target in conf.sample_targets:
        bowtie = attach_back(workflow,
            ShellCommand(
                "{tool} -p {param[threads]} -S -m {param[max_align]} \
                {param[genome_index]} {input[fastq]} {output[sam]} 2> {output[bowtie_summary]}",
                input={"genome_dir": os.path.dirname(conf.get_path("lib", "genome_index")),
                       "fastq": target + ".fastq"},
                output={"sam": target + ".sam",
                        "bowtie_summary": target + "_bowtie_summary.txt", },
                tool="bowtie",
                param={"threads": 4,
                       "max_align": 1,
                       "genome_index": conf.get_path("lib", "genome_index")}))

        bowtie.update(param=conf.items("bowtie"))
    __sam2bam(workflow, conf)

    ## using bowtie standard error output
    attach_back(workflow,
        PythonCommand(stat_bowtie,
            input={"bowtie_summaries": [t + "_bowtie_summary.txt" for t in conf.sample_targets],
                   "db": ChiLinQC_db,
                   "template": r_template},
            output={"json": conf.json_prefix + "_bowtie.json",
                    "R": conf.prefix + "_bowtie.R",
                    "pdf": conf.prefix + "_bowtie.pdf"},
            param={"sams": [t + ".sam" for t in conf.sample_targets], }))


def _bowtie_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_bowtie,
            input={"json": conf.json_prefix + "_bowtie.json",
                   "template": latex_template},
            output={"latex": conf.latex_prefix + "_bowtie.latex"}))


def __sam2bam(workflow, conf):
# convert the files from SAM to BAM format
    for target in conf.sample_targets:
        attach_back(workflow,
            ShellCommand(
                "{tool} view -bt {input[chrom_len]} {input[sam]} -o {param[tmp_bam]} && \
                {tool} sort -m {param[max_mem]} {param[tmp_bam]} {param[output_prefix]}",
                tool="samtools",
                input={"sam": target + ".sam", "chrom_len": conf.get_path("lib", "chrom_len")},
                output={"bam": target + ".bam"},
                param={"tmp_bam": target + ".tmp.bam", "output_prefix": target,
                       "max_mem": 5000000000}, )) # Use 5G memory as default
        workflow.update(param=conf.items("sam2bam"))


def _macs2(workflow, conf):
    # merge all treatments into one
    merge_bams_treat = ShellCommand(
        "{tool} merge {output[merged]} {param[bams]}",
        tool="samtools",
        input=[target + ".bam" for target in conf.treatment_targets],
        output={"merged": conf.prefix + "_treatment.bam"})
    merge_bams_treat.param = {"bams": " ".join(merge_bams_treat.input)}

    if len(conf.treatment_targets) > 1:
        attach_back(workflow, merge_bams_treat)
    else:
        # when there's only one treatment sample, use copying instead of merging
        attach_back(workflow, make_copy_command(merge_bams_treat.input[0], merge_bams_treat.output["merged"]))

    # merging step will be skipped if control sample does not exist
    # So be careful to check whether there are control samples before using `_control.bam`
    if len(conf.control_targets) > 1:
        merge_bams_control = merge_bams_treat.clone
        merge_bams_control.input = [target + ".bam" for target in conf.control_targets]
        merge_bams_control.output = {"merged": conf.prefix + "_control.bam"}
        merge_bams_control.param = {"bams": " ".join(merge_bams_control.input)}
        attach_back(workflow, merge_bams_control)
    elif len(conf.control_targets) == 1:
        attach_back(workflow, make_copy_command(conf.control_targets[0] + ".bam", conf.prefix + "_control.bam"))

    macs2_on_merged = attach_back(workflow, ShellCommand(
        "{tool} callpeak -B -q 0.01 --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel \
        {param[treat_opt]} {param[control_opt]} -n {param[description]}",
        tool="macs2",
        input={"treat": conf.prefix + "_treatment.bam"},
        output={"peaks": conf.prefix + "_peaks.bed",
                "summit": conf.prefix + "_summits.bed",
                "treat_bdg": conf.prefix + "_treat_pileup.bdg",
                "ENCODE": conf.prefix + "_peaks.encodePeak",
                "peaks_xls": conf.prefix + "_peaks.xls",
                "control_bdg": conf.prefix + "_control_lambda.bdg"},
        param={"description": conf.prefix,
               "keep_dup": 1,
               "shiftsize": 73},
        name="macs2_callpeak_merged"))
    macs2_on_merged.param["treat_opt"] = "-t " + macs2_on_merged.input["treat"]

    # control option is skipped if control samples does not exist
    if len(conf.control_targets) >= 1:
        macs2_on_merged.input["control"] = conf.prefix + "_control.bam"
        macs2_on_merged.param["control_opt"] = "-c " + macs2_on_merged.input["control"]
    else:
        macs2_on_merged.param["control_opt"] = ""

    macs2_on_merged.update(param=conf.items("macs2"))


    # For bedGraphToBigwiggle bugs, we need to remove coordinates over-border coordinates
    # As _control_lambda.bdg always exist. There are no need to check whether there are control samples.
    bdg_trim_control = attach_back(workflow,
        ShellCommand(
            '{tool} intersect -a {input[bdg]} -b {input[chrom_bed]} -wa -f 1.00 > {output}',
            tool="bedtools",
            input={"bdg": conf.prefix + "_control_lambda.bdg",
                   'chrom_bed': conf.get_path("lib", "chrom_bed")},
            output=conf.prefix + "_control_lambda.bdg.tmp",
            name="bedGraph filtering"))

    bdg_trim_treat = bdg_trim_control.clone
    bdg_trim_treat.input["bdg"] = conf.prefix + "_treat_pileup.bdg"
    bdg_trim_treat.output = conf.prefix + "_treat_pileup.bdg.tmp"
    attach_back(workflow, bdg_trim_treat)

    bdg2bw_treat = attach_back(workflow,
        ShellCommand(
            "{tool} {input[bdg]} {input[chrom_len]} {output[bw]}",
            tool="bedGraphToBigWig",
            input={"bdg": conf.prefix + "_control_lambda.bdg.tmp",
                   "chrom_len": conf.get("lib", "chrom_len")},
            output={"bw": conf.prefix + "_control.bw"},
            name="bdg_to_bw"))

    # prototype used here to do the similar thing on treatment file
    bdg2bw_control = bdg2bw_treat.clone
    bdg2bw_control.input["bdg"] = conf.prefix + "_treat_pileup.bdg.tmp"
    bdg2bw_control.output["bw"] = conf.prefix + "_treat.bw"
    attach_back(workflow, bdg2bw_control)

    attach_back(workflow, PythonCommand(
        stat_macs2,
        input={"macs2_peaks_xls": conf.prefix + "_peaks.xls",
               "db": ChiLinQC_db,
               "template": r_template},
        output={"json": conf.json_prefix + "_macs2.json",
                "R": conf.prefix + "_macs2.R",
                "pdf": conf.prefix + "_macs2.pdf"},
        param={"id": conf.id},
        name="MACS2 summary"))


def _macs2_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_macs2,
            input={"json": conf.json_prefix + "_macs2.json",
                   "template": latex_template},
            output={"latex": conf.latex_prefix + "_macs2.latex"}))


def _macs2_on_sample(workflow, conf):
    # Though macs command already exists, I choose not to use prototype here
    # Because the prototype definition and usage might be far from each other, making codes not readable

    # TODO: move the condition of rep here
    for target in conf.treatment_targets:
        macs2_on_rep = attach_back(workflow,
            ShellCommand(
                "{tool} callpeak -B -q 0.01 --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel \
                {param[treat_opt]} {param[control_opt]} -n {param[description]}",
                tool="macs2",
                input={"treat": target + ".bam"},
                output={"peaks": target + "_peaks.bed",
                        "summit": target + "_summits.bed",
                        "treat_bdg": target + "_treat_pileup.bdg",
                        "ENCODE": target + "_peaks.encodePeak",
                        "peaks_xls": target + "_peaks.xls",
                        "control_bdg": target + "_control_lambda.bdg"},
                param={"description": target, "keep_dup": 1, "shiftsize": 73},
                name="macs2_callpeak_rep"))
        macs2_on_rep.param["treat_opt"] = "-t " + macs2_on_rep.input["treat"]
        # control option is skipped if control samples does not exist
        if len(conf.control_targets) >= 1:
            macs2_on_rep.input["control"] = conf.prefix + "_control.bam"
            macs2_on_rep.param["control_opt"] = "-c " + macs2_on_rep.input["control"]
        else:
            macs2_on_rep.param["control_opt"] = ""

        macs2_on_rep.update(param=conf.items("macs2"))

        ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
        ## filter bdg file to remove over-border coordinates
        bdg_trim_controlrep = attach_back(workflow,
            ShellCommand(
                '{tool} intersect -a {input} -b {param[chrom_bed]} -wa -f 1.00 > {output}',
                tool="bedtools",
                input=target + "_control_lambda.bdg",
                output=target + "_control_lambda.bdg.tmp",
                param={'chrom_bed': conf.get("lib", "chrom_bed")},
                name="bedGraph replicate filtering"))
        bdg_trim_treatrep = bdg_trim_controlrep.clone
        bdg_trim_treatrep.input = target + "_treat_pileup.bdg"
        bdg_trim_treatrep.output = target + "_treat_pileup.bdg.tmp"
        attach_back(workflow, bdg_trim_treatrep)

        bdg2bw_treatrep = attach_back(workflow,
            ShellCommand(
                "{tool} {input} {param[chrom_len]} {output}",
                tool="bedGraphToBigWig",
                input=target + "_treat_pileup.bdg.tmp",
                output=target + "_treat.bw",
                param={"chrom_len": conf.get("lib", "chrom_len")}, name="bdg_to_bw"))

        # prototype used here to do the similar thing on treatment file
        bdg2bw_controlrep = bdg2bw_treatrep.clone
        bdg2bw_controlrep.input = target + "_treat_pileup.bdg.tmp"
        bdg2bw_controlrep.output = target + "_treat.bw"
        attach_back(workflow, bdg2bw_controlrep)

    attach_back(workflow,
        PythonCommand(
            stat_macs2_on_sample,
            input={"all_peak_xls": [target + "_peaks.xls" for target in conf.treatment_targets],
                   "db": ChiLinQC_db,
                   "template": r_template},
            output={"R": conf.prefix + "_macs2_on_sample.R",
                    "json": conf.json_prefix + "_macs2_on_sample.json",
                    "pdf": conf.prefix + "_macs2_on_sample.pdf"},
            param={"ids": conf.treatment_bases}))


def _macs2_on_sample_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_macs2_on_sample,
            input={"json": conf.json_prefix + "_macs2_on_sample.json",
                   "template": latex_template},
            output={"latex": conf.latex_prefix + "_macs2_on_sample.latex"}))


def _macs2_venn(workflow, conf):
    # awk and bedClip to remove outlier for venn and correlation plot
    for target in conf.treatment_targets:
        bed_filter = attach_back(workflow,
            ShellCommand(
                "{tool} '{{if ($2 >= 0 && $2 < $3) print}}' {input} > {output}",
                tool="awk",
                input=target + "_peaks.bed",
                output=target + "_peaks.bed.tmp",
                name="filter bed files"))

        # prototype used here to do the similar thing on bedclip
        bed_clip = attach_back(workflow,
            ShellCommand(
                template="{tool} {input} {param[chrom_len]} {output}",
                tool="bedClip",
                input=target + "_peaks.bed.tmp",
                output=target + "_peaks.bed",
                param={'chrom_len': conf.get_path("lib", "chrom_len")},
                name="bedclip filter"))
        bed_clip.allow_fail = True

    venn_on_peaks = attach_back(workflow,
        ShellCommand(
            "{tool} -t Overlap_of_Replicates {param[beds]} && \
            mv venn_diagram.png {output}",
            tool="venn_diagram.py",
            input=[target + "_peaks.bed.tmp" for target in conf.treatment_targets],
            output=conf.prefix + "_venn.png", name="venn_diagram"))
    venn_on_peaks.param = {"beds": " ".join(venn_on_peaks.input)}
    venn_on_peaks.allow_fail = True

    attach_back(workflow,
        PythonCommand(
            stat_venn,
            input={"venn": conf.prefix + "_venn.png"},
            output={"json": conf.json_prefix + "_venn.json"}))


def _macs2_venn_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_venn,
            input={"json": conf.json_prefix + "_venn.json",
                   "template": latex_template},
            output={"latex": conf.latex_prefix + "_venn.latex"}))


def _macs2_cor(workflow, conf):
    cor_on_bw = attach_back(workflow,
        ShellCommand(
            template=
            """{tool} \
            -s {param[wig_correlation_step]}  \
            --min-score {param[wig_correlation_min]} --max-score {param[wig_correlation_max]} \
            -r {output[R]} {param[bw]}  {param[rep]} && \
            mv {output[R]}.pdf {output[pdf]}""",
            tool="bigwig_correlation.py",
            input=[target + "_treat.bw" for target in conf.treatment_targets],
            output={"R": conf.prefix + "_cor.R", "pdf": conf.prefix + "_cor.pdf"},
            param={"wig_correlation_method": "mean",
                   "wig_correlation_min": 2,
                   "wig_correlation_max": 50},
            name="cor_on_bw"))
    cor_on_bw.param["bw"] = " ".join(cor_on_bw.input)
    cor_on_bw.param["rep"] = " ".join([" -l replicate_%s" % (x + 1) for x in range(len(conf.treatment_pairs))])
    cor_on_bw.update(param=conf.items("correlation"))
    cor_on_bw.allow_fail = True

    attach_back(workflow,
        PythonCommand(
            stat_cor,
            input={"correlation_R": conf.prefix + "_cor.R",
                   "cor_pdf": conf.prefix + "_cor.pdf"},
            output={"json": conf.json_prefix + "_cor.json"}))


def _macs2_cor_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_cor,
            input={"json": conf.json_prefix + "_cor.json",
                   "template": latex_template},
            output={"latex": conf.latex_prefix + "_cor.latex"}))


def _DHS_overlap(workflow, conf):
    attach_back(workflow,
        ShellCommand(
            "{tool} -wa -u  \
            -a {input[MACS2_bed]} -b {input[DHS_peaks_bed]} > {output}",
            tool="intersectBed",
            input={"MACS2_bed": conf.prefix + "_peaks.bed",
                   "DHS_peaks_bed": conf.get("lib", "dhs")},
            output=conf.prefix + "_DHS_overlap_peaks_bed",
            param=None))

    attach_back(workflow, PythonCommand(
        stat_dhs,
        input={"dhs_peaks": conf.prefix + "_DHS_overlap_peaks_bed",
               "macs2_peaks_xls": conf.prefix + "_peaks.xls"},
        output={"json": conf.json_prefix + "_dhs.json"},
        param=None,
        name="DHS summary"))


def _velcro_overlap(workflow, conf):
    attach_back(workflow,
        ShellCommand(
            "{tool} -wa -u  \
            -a {input[MACS2_bed]} -b {input[velcro_peaks_bed]} > {output}",
            tool="intersectBed",
            input={"MACS2_bed": conf.prefix + "_peaks.bed",
                   "velcro_peaks_bed": conf.get("lib", "velcro")},
            output=conf.prefix + "_velcro_overlap_peaks_bed",
            param=None))

    attach_back(workflow, PythonCommand(
        stat_velcro,
        input={"velcro_peaks": conf.prefix + "_velcro_overlap_peaks_bed",
               "macs2_peaks_xls": conf.prefix + "_peaks.xls"},
        output={"json": conf.json_prefix + "_velcro.json"},
        param=None,
        name="Velcro summary"))


def _ceas(workflow, conf):
    get_top_peaks = attach_back(workflow,
        ShellCommand(
            "{tool} -r -g -k 5 {input} | head -n {param[peaks]} > {output}",
            tool="sort",
            input=conf.prefix + "_peaks.bed",
            output=conf.prefix + "_peaks_top.bed",
            param={"peaks": 3000},
            name="top summits for conservation"))
    get_top_peaks.update(param=conf.items("ceas"))

    ceas = attach_back(workflow,
        ShellCommand(
            "{tool} -g {input[gene_table]} \
            -l {input[chrom_len]} \
            --name {param[description]} -b {input[bed]} -w {input[bw]}",
            tool="ceasBW",
            input={"bed": conf.prefix + "_peaks.bed",
                   "bw": conf.prefix + "_treat.bw",
                   "gene_table": conf.get_path("lib", "gene_table"),
                   "chrom_len": conf.get_path("lib", "chrom_len")},
            output={"R": conf.prefix + "_ceas.R",
                    "pdf": conf.prefix + "_ceas.pdf"},
            param={"description": conf.prefix + "_ceas"},
            # `description` is the prefix of output, so "_ceas" here is need
            name="ceas"))
    ceas.update(param=conf.items("ceas"))

    attach_back(workflow,
        PythonCommand(
            stat_ceas,
            input={"macs2_peaks_xls": conf.prefix + "_peaks.xls",
                   "ceas_rscript": conf.prefix + "_ceas.R", },
            output={"R": conf.prefix + "_qc_ceas.R",
                    "peakheight_and_pie_pdf": conf.prefix + "_peakheight_and_pie.pdf",
                    "metagene_dist_pdf": conf.prefix + "_metagene_dist.pdf",
                    "json": conf.json_prefix + "_ceas.json"}, ))


def _ceas_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_ceas,
            input={"json": conf.json_prefix + "_ceas.json",
                   "template": latex_template},
            output={"latex": conf.latex_prefix + "_ceas.latex"}))


def _conservation(workflow, conf):
    get_top_peaks = attach_back(workflow,
        ShellCommand(
            "{tool} -r -g -k 5 {input} | head -n {param[peaks]} > {output}",
            tool="sort",
            input=conf.prefix + "_summits.bed",
            output=conf.prefix + "_summits_topconserv.bed",
            param={'peaks': 5000}, name="top summits for conservation"))
    get_top_peaks.update(param=conf.items('conservation'))

    conservation = attach_back(workflow,
        ShellCommand(
            "{tool} -t Conservation_at_summits \
            -d {input[phast]} -l Peak_summits {input[bed]} -w {param[width]} &&\
            mv tmp.pdf {output[pdf]} && \
            mv tmp.R {output[R]}",

            tool="conservation_plot.py",
            input={"bed": conf.prefix + "_summits_topconserv.bed",
                   "phast": conf.get_path("lib", "phast")},
            output={"pdf": conf.prefix + "_conserv_tmp.pdf",
                    "R": conf.prefix + "_conserv_tmp.R"},
            param={"width": 4000},
            name="conservation"))
    conservation.update(param=conf.items('conservation'))

    attach_back(workflow,
        ShellCommand(
            "{tool} -resize 500x500 -density 50  {input[pdf]} {output[pdf]} && mv {input[R]} {output[R]}",
            tool="convert", ## width differs histone mark and TF
            input={"pdf": conf.prefix + "_conserv_tmp.pdf",
                   "R": conf.prefix + "_conserv_tmp.R"},
            output={"pdf": conf.prefix + "conserv.pdf", "R": conf.prefix + "conserv.R"},
            name="convert pdf to png", ))

    attach_back(workflow,
        PythonCommand(
            stat_conservation,
            input={"conservationR": conf.prefix + "conserv.R",
                   "historical_conservation_cluster_text": resource_filename("chilin2", "db/TF_centers.txt")},
            output={"R": conf.prefix + "_qc_conserv_compare.R",
                    "compare_pdf": conf.prefix + "_qc_conserv_compare.pdf",
                    "pdf": conf.prefix + "conserv.pdf",
                    "json": conf.json_prefix + "_conserv.json"},
            param={"atype": conf.get("Basis", "type", "TF"), "id": conf.id}))


def _conservation_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_conservation,
            input={"json": conf.json_prefix + "_conserv.json",
                   "template": latex_template},
            output={"latex": conf.latex_prefix + "_conserv.latex"}))


def _seqpos(workflow, conf):
    # This will work for only Human and Mouse
    # MDseqpos take uniform input, we need to remove random sequence part and get top 1000 summits
    get_top_summits = attach_back(workflow,
        ShellCommand(
            '{tool} "/^chr[1-22XY]/" {input} |sort -r -g -k 5|head -n {param[peaks]} > {output}',
            tool="awk",
            input=conf.prefix + "_summits.bed",
            output=conf.prefix + "_summits_topmdfilter.bed",
            param={'peaks': 1000}, name="filter summits for motif"))
    get_top_summits.update(param=conf.items("seqpos"))

    mdseqpos = attach_back(workflow,
        ShellCommand(
            "{tool} -d  -w 600  -p 0.001  -m cistrome.xml  -s {param[species]} {input} {param[version]}",
            "MDSeqPos.py",
            input=conf.prefix + "_summits_topmdfilter.bed",
            output="results",
            param={"species": "hs", "version": "hg19"}, name="motif finding"))
    mdseqpos.update(param=conf.items("seqpos"))

    attach_back(workflow,
        ShellCommand(
            "{tool} {input} {output}",
            "mv",
            input="results",
            output=conf.prefix + "_seqpos",
            name="mv seqpos"))

    attach_back(workflow,
        PythonCommand(
            stat_seqpos,
            input={"seqpos": conf.prefix + "_seqpos/" + "mdseqpos_out.html"},
            output={"json": conf.json_prefix + "_seqpos.json"},
            param={"z_score_cutoff": -15}))

def prepare_summary_plain(workflow, conf):
    """summary of all criteria by plain text"""

    return

def _seqpos_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_seqpos,
            input={"json": conf.json_prefix + "_seqpos.json",
                   "template": latex_template},
            output={"latex": conf.latex_prefix + "_seqpos.latex"}))


def _ending_latex(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            latex_end,
            input={"template": latex_template},
            output={"latex": conf.latex_prefix + "_ending.latex"}))

def prepare_clean_up(workflow, conf):
    """
    package all the necessary results and delete temporary files
    """
    bams = glob('*.bam')
    xls = glob('*.xls')
    summits = glob('*_summits.bed')
    peaks = glob('*_peaks.bed')
    bw = glob('*.bw')
    png = glob('*.png')
    cor = glob('*cor*')
    pdf = glob('*_ceas.pdf')
    r = glob('*_ceas.R')
    m = glob('*.zip')
    su = glob('dataset*.txt')
    qc = glob('*.tex')
    qcp = glob('*QC.pdf')
    fls = [bams, xls, summits, peaks, bw, png, pdf, r, m, cor, su, qc, qcp]
    folder = conf.target_dir + '/dataset' + conf.id
    os.system('mkdir %s' % folder)
    for fs in fls:
        os.system('cp %s %s' % (' '.join(fs), folder))
    preservation = os.listdir(".")
    for f in preservation:
        if f.startswith("dataset") and os.path.isdir(f):
            continue
        else:
            os.system("rm -r %s" % f)

class StepChecker:
    def __init__(self, start, end, skips):
        self.start = start
        self.end = end
        self.skips = skips

    def need_run(self, step_id):
        if step_id < self.start:
            return False
        if step_id > self.end:
            return False
        if step_id in self.skips:
            return False
        return True


class ChiLinBuilder:
    def __init__(self, workflow, conf):
        self.workflow = workflow
        self.conf = conf
        self.LaTex_fragments = []
        self.plain_fragments = []
        self.finished = set()

    def build(self, prep_func, tag=None):
        prep_func(self.workflow, self.conf)
        if tag:
            self.finished.add(tag)

    def attach_back(self, command):
        attach_back(self.workflow, command)


def create_workflow(args, conf, step_checker : StepChecker):
    """
    :type conf:ChiLinConfig
    """



    # Whether there are replicates for treatment group
    have_treat_reps = len(conf.treatment_pairs) >= 2
    has_dhs = conf.get("lib", "dhs")
    has_velcro = conf.get("lib", "velcro")
    need_run = step_checker.need_run

    bld = ChiLinBuilder(Workflow(name="Main"), conf)

    bld.attach_back(ShellCommand(
        "if [ ! -d '{output}' ]; then mkdir -p {output}; fi",
        output=conf.target_dir))

    if need_run(1):
        bld.build(_groom_sequencing_files)

    if need_run(2):
        bld.build(_raw_QC)
        bld.build(_raw_QC_latex)
        #TODO: merge
        bld.build(prepare_library_contamination)
        bld.build_LaTex(prepare_library_contaminationTex)

    if need_run(3):
        bld.build(_bowtie)
        bld.build(_bowtie_latex)

    if need_run(4):
        bld.build(_macs2)
        bld.build(_macs2_latex)

    if need_run(5):
        bld.build(_macs2_on_sample)
        bld.build(_macs2_on_sample_latex)

    if have_treat_reps:
        if need_run(6):
            bld.build(_macs2_venn)
            bld.build(_macs2_venn_latex)

        if need_run(7):
            bld.build(_macs2_cor)
            bld.build(_macs2_cor_latex)

    if has_dhs and need_run(8):
        bld.build(_DHS_overlap)

    if has_velcro and need_run(9):
        bld.build(_velcro_overlap)

    if need_run(10):
        bld.build(_ceas)
        bld.build(_ceas_latex)

    if need_run(11):
        bld.build(_conservation)
        bld.build(_conservation_latex)

    if need_run(12):
        bld.build(_seqpos)
        bld.build(_seqpos_latex)

    bld.build(_ending_latex)

    if args.clean:
        bld.build(prepare_clean_up)
    return bld.workflow


def main(args=None):
    args = parse_args(args)
    print("Arguments:", args)

    conf = ChiLinConfig(args.config)
    if args.skip_step:
        skipped_steps = [int(i) for i in args.skip_step.split(",")]
    else:
        skipped_steps = []

    step_checker = StepChecker(args.start_step, args.end_step, skipped_steps)

    workflow = create_workflow(args, conf, step_checker)

    workflow.set_option(
        verbose_level=args.verbose_level,
        dry_run_mode=args.dry_run,
        resume=args.resume,
        allow_dangling=args.allow_dangling)

    workflow.invoke()

if __name__ == "__main__":
    main()
