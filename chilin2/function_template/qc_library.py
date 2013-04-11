import random
import sys
from chilin2.function_template.qc_bowtie import _qc_bowtie_summary_parse, underline_to_space
from chilin2.jinja_template_render import JinjaTemplateCommand, write_into

__author__ = 'qinqianhappy'

def write_random_records(input = {"fastq": ""}, output = {"fastq_sample": ""}, param = {"random_number": ""}):
    """ get N random headers from a fastq file without reading the
    whole thing into memory"""
    records = sum(1 for _ in open(input["fastq"])) / 4
    rand_records = sorted([random.randint(0, records - 1) for _ in range(int(param["random_number"]))])
    fha = open(input["fastq"],"rU")
    suba = open(output["fastq_sample"], "w")

    rec_no = - 1
    for rr in rand_records:
        while rec_no < rr:
            rec_no += 1
            for i in range(4): fha.readline()
        for i in range(4):
            suba.write(fha.readline())
            rec_no += 1
    suba.close()
    fha.close()

    print("wrote to %s" % (output["fastq_sample"]))

## summary of library contamination
def library_summary(input = {"latex_template": "", "bowtie_summary": "", "other1_bowtie_summary": "", "other2_bowtie_summary": ""},
                    output = {"latex_section": ""}, param = {"samples": "","species": "", "other1": "", "other2": ""}):

    bowtie_summaries = {param["species"]: [],
                        param["other1"]: [],
                        param["other2"]: []}

    for a_summary in input["bowtie_summary"]:
        parsed_summary = _qc_bowtie_summary_parse(a_summary)
        bowtie_summaries[param["species"]].append(parsed_summary["mappable_rate"])
    for a_summary in input["other1_summary"]:
        parsed_summary = _qc_bowtie_summary_parse(a_summary)
        bowtie_summaries[param["other1"]].append(parsed_summary["mappable_rate"])
    for a_summary in input["other2_summary"]:
        parsed_summary = _qc_bowtie_summary_parse(a_summary)
        bowtie_summaries[param["other2"]].append(parsed_summary["mappable_rate"])
        # basic mappable table, two layer list

    # COL 1 samples
    # COL 2 human mappable ratio (suppose correct species is human)
    # COL 3 mouse mappable ratio
    # COL 4 rat mappable ratio
    library_contamination = [["sample_name", param["species"], param["other1"], param["other2"]]]
    library_contamination = zip(map(underline_to_space, param["samples"]), bowtie_summaries[param["species"]], bowtie_summaries[param["other1"]], bowtie_summaries[param["other2"]])
    library_quality_latex = JinjaTemplateCommand(
        name="library contamination",
        template=input["latex_template"],
        param={"section_name": "library_contamination",
               "qc_report_begin": True,
               "library_contamination": library_contamination,
               'prefix_dataset_id': param['id']})
    write_into(library_quality_latex, output['latex_section'])
