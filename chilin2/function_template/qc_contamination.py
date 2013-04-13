import random
import sys
from chilin2.function_template.qc_bowtie import _qc_bowtie_summary_parse, underline_to_space
from chilin2.helpers import JinjaTemplateCommand, template_dump

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
def library_summary(input = {"latex_template": ""},
                    output = {"latex_section": ""}, param = {"samples": "", "bowtie_summaries": "", "species": "", "id": ""}):

    library_contamination = [["sample_name"] + param["species"]]

    # COL 1 samples
    # COL 2 species1 mappable ratio (suppose correct species is human)
    # COL 3 species2 mappable ratio
    # COL 4 species3 mappable ratio
    for a_summary, s in zip(param["bowtie_summary"], map(underline_to_space, param["samples"])):
        ## each bowtie_summary has several species information
        bowtie_summaries = []


        for i in a_summary:
            print(i)
            ## species 1, species2, species3
            parsed_summary = _qc_bowtie_summary_parse(i)
            bowtie_summaries.append(s)
            bowtie_summaries.append(parsed_summary["mappable_rate"])
        library_contamination.append(bowtie_summaries)

    library_quality_latex = JinjaTemplateCommand(
        name="library contamination",
        template=input["latex_template"],
        param={"section_name": "library_contamination",
               "qc_report_begin": True,
               "library_contamination": library_contamination,
               'prefix_dataset_id': param['id'],
               "render_dump": output["latex_section"]
               })
    template_dump(library_quality_latex)
