#!/usr/bin/env python3

import os
import sys
import re
import argparse
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import Workflow, attach_back
from chilin2.config import ChiLinConfig
from pkg_resources import resource_filename
from chilin2.function_template.qc_bowtie import qc_bowtie_summary_draw
from chilin2.function_template.qc_fastqc_raw_sequence import python_fastqc_dist_draw
from chilin2.function_template.qc_macs2 import qc_high_confident_peaks_draw, qc_non_redundant_rate_draw
from chilin2.function_template.qc_venn_replicate import qc_replicate_parse, qc_venn
from chilin2.function_template.qc_ceas import qc_redraw_ceas_graph
from chilin2.function_template.qc_phast_conservation import qc_conservation_draw
from chilin2.function_template.qc_mdseqpos import qc_mdseqpos_parse_and_filter_by_z_score, end_tex
from chilin2.function_template.text_summary import text_bowtie_summary, macs2_summary, dhs_summary, velcro_summary

ChiLinQC_db = resource_filename("chilin2", "db/ChiLinQC.db")
R_cumulative_template = resource_filename("chilin2", "jinja_template/R_culmulative_plot.R.jinja2")
Latex_summary_report_template = resource_filename("chilin2", "jinja_template/Latex_summary_report.jinja2")
text_summary_report_template = resource_filename("chilin2", "jinja_template/text_summary_report.jinja2")


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


def prepare_groom_sequencing_files(workflow, conf):
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
    return not_groomed


def prepare_raw_QC(workflow, conf):
    for target in conf.sample_targets:
        fastqc_run = attach_back(workflow,
            ShellCommand(
                "{tool} {input} --extract -t {param[threads]} -o {output[target_dir]}",
                input=target + ".fastq",
                output={"target_dir": conf.target_dir,
                        "summary": target + "_fastqc/fastqc_data.txt"},
                tool="fastqc",
                param={"threads": 4}))
        fastqc_run.update(param=conf.items("fastqc"))


def prepare_raw_QC_Tex(workflow, conf):
    Tex_step = attach_back(workflow,
        PythonCommand(
            python_fastqc_dist_draw,
            input={"db": ChiLinQC_db,
                   "fastqc_summary_list": [target + "_fastqc/fastqc_data.txt" for target in conf.sample_targets],
                   "R_template": R_cumulative_template,
                   "latex_template": Latex_summary_report_template},
            output={"rfile": conf.prefix + "_raw_sequence_qc.R",
                    "latex_section": conf.prefix + "_raw_sequence_qc.tex",
                    "pdf": conf.prefix + "_raw_sequence_qc.pdf"},
            param={"ids": conf.sample_bases,
                   "id": conf.id}))
    return Tex_step.output["latex_section"]


def prepare_bowtie(workflow, conf):
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
    _prepare_sam2bam(workflow, conf)


def prepare_bowtie_plain(workflow, conf):
    ## using bowtie standard error output 
    attach_back(workflow,
        PythonCommand(text_bowtie_summary,
            input={"data_template": text_summary_report_template,
                   "bowtie_summary": [t + "_bowtie_summary.txt" for t in conf.sample_targets]},
            output={"sum_section": conf.prefix + "_all_bowtie_summary"},
            param={"sam_files": [t + ".sam" for t in conf.sample_targets], }))
    return conf.prefix + "_all_bowtie_summary"


def prepare_bowtie_Tex(workflow, conf):
    Tex_step = attach_back(workflow, PythonCommand(qc_bowtie_summary_draw,
        input={"all_bowtie_summary": [target + "_bowtie_summary.txt" for target in conf.sample_targets],
               "R_template": R_cumulative_template,
               "latex_template": Latex_summary_report_template,
               "db": ChiLinQC_db},
        output={"rfile": conf.prefix + "pable_rate.R",
                "pdf": conf.prefix + "pable_rate.pdf",
                "latex_section": conf.prefix + "pable.tex"},
        param={"ids": conf.sample_bases}))
    return Tex_step.output["latex_section"]


def _prepare_sam2bam(workflow, conf):
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


def prepare_macs2(workflow, conf):
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


def prepare_macs2_Tex(workflow, conf):
    Tex_step = attach_back(workflow, PythonCommand(
        qc_high_confident_peaks_draw,
        input={"macs2_peaks_xls": conf.prefix + "_peaks.xls",
               "R_template": R_cumulative_template,
               "latex_template": Latex_summary_report_template,
               "db": ChiLinQC_db},
        output={"rfile": conf.prefix + "_high_confident_peak_rate.R",
                "latex_section": conf.prefix + "_high_confident.tex",
                "pdf": conf.prefix + "_high_confident_peak_rate.pdf"},
        param={"id": os.path.basename(conf.prefix)},
        name="high_confident_peaks"
    ))
    return Tex_step.output["latex_section"]


def prepare_macs2_on_rep(workflow, conf):
    # Though macs command already exists, I choose not to use prototype here
    # Because the prototype definition and usage might be far from each other, making codes not readable

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


def prepare_macs2_on_rep_Tex(workflow, conf):
    Tex_step = attach_back(workflow,
        PythonCommand(
            qc_non_redundant_rate_draw,
            input={"all_peak_xls": [target + "_peaks.xls" for target in conf.treatment_targets],
                   "db": ChiLinQC_db,
                   "R_template": R_cumulative_template,
                   "latex_template": Latex_summary_report_template},
            output={"rfile": conf.prefix + "_redundant_dist.R",
                    "latex_section": conf.prefix + "_redundant.tex",
                    "pdf": conf.prefix + "_redundant_dist.pdf"},
            param={"ids": conf.treatment_bases}))
    return Tex_step.output["latex_section"]


def prepare_macs2_venn_on_rep(workflow, conf):
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


def prepare_macs2_venn_on_rep_Tex(workflow, conf):
    Tex_step = attach_back(workflow,
        PythonCommand(
            qc_venn,
            input={"venn": conf.prefix + "_venn.png",
                   "latex_template": Latex_summary_report_template},
            output={"latex_section": conf.prefix + "_venn.tex"}))
    return Tex_step.output["latex_section"]


def prepare_macs2_cor_on_rep(workflow, conf):
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


def prepare_macs2_cor_on_rep_Tex(workflow, conf):
    Tex_step = attach_back(workflow,
        PythonCommand(
            qc_replicate_parse,
            input={"correlation_R": conf.prefix + "_cor.R",
                   "latex_template": Latex_summary_report_template,
                   "cor_pdf": conf.prefix + "_cor.pdf"},
            output={"latex_section": conf.prefix + "_cor.tex"}))
    return Tex_step.output["latex_section"]


def prepare_DHS_overlap(workflow, conf):
    DHS_overlap = attach_back(workflow,
        ShellCommand(
            "{tool} -wa -u  \
            -a {input[MACS2_bed]} -b {input[DHS_peaks_bed]} > {output}",
            tool="intersectBed",
            input={"MACS2_bed": conf.prefix + "_peaks.bed",
                   "DHS_peaks_bed": conf.get("lib", "dhs")},
            output=conf.prefix + "_DHS_overlap_peaks_bed",
            param=None))


def prepare_velcro_overlap(workflow, conf):
    velcro_overlap = attach_back(workflow,
        ShellCommand(
            "{tool} -wa -u  \
            -a {input[MACS2_bed]} -b {input[velcro_peaks_bed]} > {output}",
            tool="intersectBed",
            input={"MACS2_bed": conf.prefix + "_peaks.bed",
                   "velcro_peaks_bed": conf.get("lib", "velcro")},
            output=conf.prefix + "_velcro_overlap_peaks_bed",
            param=None))


def prepare_macs2_plain(workflow, conf):
    # using macs2 merged peaks.xls information for unique location and peaks summary
    # DHS peaks and Velcro peaks summary
    plain_step = attach_back(workflow, PythonCommand(
        macs2_summary,
        input={"data_template": text_summary_report_template,
               "macs2_peaks_xls": conf.prefix + "_peaks.xls"},
        output={"sum_section": conf.prefix + "_peaks_summary"},
        param=None,
        name="MACS2 summary"))
    return plain_step.output["sum_section"]


def prepare_peaks_dhs_plain(workflow, conf):
    plain_step = attach_back(workflow, PythonCommand(
        dhs_summary,
        input={"data_template": text_summary_report_template, "dhs_peaks": conf.prefix + "_DHS_overlap_peaks_bed",
               "macs2_peaks_xls": conf.prefix + "_peaks.xls"},
        output={"sum_section": conf.prefix + "_dhs_summary"},
        param=None,
        name="DHS summary"))
    return plain_step.output["sum_section"]


def prepare_peaks_velcro_plain(workflow, conf):
    plain_step = attach_back(workflow, PythonCommand(
        velcro_summary,
        input={"data_template": text_summary_report_template, "velcro_peaks": conf.prefix + "_velcro_overlap_peaks_bed",
               "macs2_peaks_xls": conf.prefix + "_peaks.xls"},
        output={"sum_section": conf.prefix + "_velcro_summary"},
        param=None,
        name="Velcro summary"))
    return plain_step.output["sum_section"]


def prepare_ceas(workflow, conf):
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


def prepare_ceas_Tex(workflow, conf):
    Tex_step = attach_back(workflow,
        PythonCommand(
            qc_redraw_ceas_graph,
            input={"macs2_peaks_xls": conf.prefix + "_peaks.xls",
                   "ceas_rscript": conf.prefix + "_ceas.R",
                   "latex_template": Latex_summary_report_template},
            output={"rfile": conf.prefix + "_qc_ceas.R",
                    "peakheight_and_pie_pdf": conf.prefix + "_peakheight_and_pie.pdf",
                    "metagene_dist_pdf": conf.prefix + "_metagene_dist.pdf",
                    "latex_section": conf.prefix + "_ceas_qc.tex",
            }, ))
    return Tex_step.output["latex_section"]


def prepare_phast_conservation(workflow, conf):
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


def prepare_phast_conservation_Tex(workflow, conf):
    Tex_step = attach_back(workflow,
        PythonCommand(
            qc_conservation_draw,
            input={"conservationR": conf.prefix + "conserv.R",
                   "historical_conservation_cluster_text": resource_filename("chilin2", "db/TF_centers.txt"),
                   "latex_template": Latex_summary_report_template},
            output={"rfile": conf.prefix + "_qc_conserv_compare.R",
                    "compare_pdf": conf.prefix + "_qc_conserv_compare.pdf",
                    "pdf": conf.prefix + "conserv.pdf",
                    "latex_section": conf.prefix + "_conserv_qc.tex"},
            param={"atype": conf.get("Basis", "type", "TF"), "id": conf.id})
    )

    return Tex_step.output["latex_section"]


def prepare_mdseqpos(workflow, conf):
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


def prepare_mdseqpos_Tex(workflow, conf):
    Tex_step = attach_back(workflow,
        PythonCommand(
            qc_mdseqpos_parse_and_filter_by_z_score,
            input={"seqpos": conf.prefix + "_seqpos/" + "mdseqpos_out.html",
                   "latex_template": Latex_summary_report_template},
            output={"latex_section": conf.prefix + "_seqpos.tex"},
            param={"z_score_cutoff": -15}))
    return Tex_step.output["latex_section"]


def prepare_Tex_ending(workflow, conf):
    Tex_step = attach_back(workflow,
        PythonCommand(end_tex,
            input=Latex_summary_report_template,
            output={"latex_section": conf.prefix + "_end.tex"},
            param=None))
    return Tex_step.output["latex_section"]


def prepare_report_summary(workflow, conf, latex_combined, text_combined):
    cat_text = attach_back(workflow,
        ShellCommand(
            "{tool} {param[input]} > {output}",
            tool="cat",
            input=text_combined,
            output=conf.prefix + "_summary.txt",
            param={"input": ' '.join(text_combined)},
            name="cat text summary"))

    cat_tex = attach_back(workflow,
        ShellCommand(
            "{tool} {param[input]} > {output}",
            input=latex_combined,
            output=conf.prefix + "_report.tex",
            param={"input": " ".join(latex_combined)},
            tool="cat",
            name="cat latex"))

    report = attach_back(workflow,
        ShellCommand(
            # Somehow the pdflatex has to be invoked twice..
            "{tool} -output-directory {output[dir]} -jobname={param[name]} {input} \
            && {tool} -output-directory {output[dir]} -jobname={param[name]} {input}",

            tool="pdflatex",
            input=conf.prefix + "_report.tex",
            # output[pdf] should use "conf.prefix" to have the absolute path
            output={"dir": conf.target_dir, "pdf": conf.prefix + "_report.pdf"},
            # param[name] should use "conf.id" to avoid using absolute path
            param={"name": conf.id + "_report"},
            name="report"))


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

    def build(self, prep_func):
        prep_func(self.workflow, self.conf)

    def attach_back(self, command):
        attach_back(self.workflow, command)

    def build_LaTex(self, Tex_prep_func):
        self.LaTex_fragments.append(Tex_prep_func(self.workflow, self.conf))

    def build_plain(self, plain_prep_func):
        self.plain_fragments.append(plain_prep_func(self.workflow, self.conf))

    def build_summary(self, summary_prep_func):
        summary_prep_func(self.workflow, self.conf, self.LaTex_fragments, self.plain_fragments)


def create_workflow(args, conf, step_checker : StepChecker):
    """
    :type conf:ChiLinConfig
    """



    # Whether there are replicates for treatment group
    have_reps = len(conf.treatment_pairs) >= 2
    has_dhs = conf.get("lib", "dhs")
    has_velcro = conf.get("lib", "velcro")
    need_run = step_checker.need_run

    bld = ChiLinBuilder(Workflow(name="Main"), conf)

    bld.attach_back(ShellCommand(
        "if [ ! -d '{output}' ]; then mkdir -p {output}; fi",
        output=conf.target_dir))

    if need_run(1):
        bld.build(prepare_groom_sequencing_files)

    if need_run(2):
        bld.build(prepare_raw_QC)

        bld.build_LaTex(prepare_raw_QC_Tex)

    if need_run(3):
        bld.build(prepare_bowtie)

        bld.build_LaTex(prepare_bowtie_Tex)
        bld.build_plain(prepare_bowtie_plain)

    if need_run(4):
        bld.build(prepare_macs2)

        bld.build_LaTex(prepare_macs2_Tex)
        bld.build_plain(prepare_macs2_plain)

    if have_reps:
        if need_run(5):
            bld.build(prepare_macs2_on_rep)
            bld.build_LaTex(prepare_macs2_Tex)

        if need_run(6):
            bld.build(prepare_macs2_venn_on_rep)
            bld.build_LaTex(prepare_macs2_venn_on_rep_Tex)

        if need_run(7):
            bld.build(prepare_macs2_cor_on_rep)
            bld.build_LaTex(prepare_macs2_cor_on_rep_Tex)

    if has_dhs and need_run(8):
        bld.build(prepare_DHS_overlap)

        bld.build_plain(prepare_peaks_dhs_plain)

    if has_velcro and need_run(9):
        bld.build(prepare_velcro_overlap)
        bld.build_plain(prepare_peaks_velcro_plain)

    if need_run(10):
        bld.build(prepare_ceas)
        bld.build_LaTex(prepare_ceas_Tex)

    if need_run(11):
        bld.build(prepare_phast_conservation)
        bld.build_LaTex(prepare_phast_conservation_Tex)

    if need_run(12):
        bld.build(prepare_mdseqpos)
        bld.build_LaTex(prepare_mdseqpos_Tex)

    bld.build_LaTex(prepare_Tex_ending)

    bld.build_summary(prepare_report_summary)


    if args.clean:
        print("test for clean")
        #clean_up()
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
