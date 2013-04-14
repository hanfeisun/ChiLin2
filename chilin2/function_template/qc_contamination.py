import json
import random
import sys
from chilin2.function_template.qc_bowtie import _bowtie_summary_parse
from chilin2.helpers import JinjaTemplateCommand, template_dump, json_load, underline_to_space

__author__ = 'qinqianhappy'

def fastq_sampling(input = {"fastq": ""}, output = {"fastq_sample": ""}, param = {"random_number": ""}):
    """ get N random headers from a fastq file without reading the
    whole thing into memory"""
    num_lines = sum(1 for _ in open(input["fastq"])) / 4

    rand_nums = sorted([random.randint(0, num_lines - 1) for _ in range(param["random_number"])])

    fastq = open(input["fastq"],"rU")
    fastq_sample = open(output["fastq_sample"], "w")

    cur_num = - 1
    for rand_num in rand_nums:
        while cur_num < rand_num:
            for i in range(4):
                fastq.readline()
            cur_num += 1

        for i in range(4):
            fastq_sample.write(fastq.readline())
        cur_num += 1

    fastq_sample.close()
    fastq.close()

    print("wrote to %s" % (output["fastq_sample"]))

## summary of library contamination
def stat_contamination(input = {"bowtie_summaries": [[]]},
                    output = {"json": ""}, param = {"samples": "", "species": "", "id": ""}):


    library_contamination = [["sample_name"] + param["species"]]


    for a_summary, s in zip(input["bowtie_summaries"], map(underline_to_space, param["samples"])):
        ## each bowtie_summary has several species information
        bowtie_summaries = []
        for i in a_summary:
            print(i)
            ## species 1, species2, species3
            parsed_summary = _bowtie_summary_parse(i)
            bowtie_summaries.append(s)
            bowtie_summaries.append(parsed_summary["mappable_rate"])
        library_contamination.append(bowtie_summaries)

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict["stat"] = library_contamination
    with open(output["json"], "w") as f:
        json.dump(json_dict, f, indent=4)


def latex_contamination(input, output, param):

    json_dict = json_load(input["json"])

    library_quality_latex = JinjaTemplateCommand(
        name="library contamination",
        template=input["template"],
        param={"section_name": "library_contamination",
               "library_contamination": json_dict["stat"],
               'prefix_dataset_id': json_dict["param"]['id'],
               "render_dump": output["latex"]
        })
    template_dump(library_quality_latex)
