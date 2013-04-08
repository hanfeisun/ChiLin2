from chilin2.jinja_template_render import JinjaTemplateCommand, write_into

def text_bowtie_summary(input={"data_template": "", "bowtie_summary": []},
                        output={"sum_section": ""},
                        param={"sam_files": []}):
    """ sams = [{'name1':a, 'total1': 5...}, {'name2':c, 'total2': 3...}...] **args """

    # unique location is in text_macs2_summary part
    sams_info = []
    for a_summary, a_sam in zip(input["bowtie_summary"], param["sam_files"]):
        summary_lines = open(a_summary).readlines()
        # first line should be like this: `# reads processed: 1234`
        total = int(summary_lines[0].strip().split(":")[1])

        # fifth line should be like this: `Reported 1234 alignments to 1 output stream(s)`
        align = int(summary_lines[4].strip().split()[1])

        usable_percentage = float(align) / total

        sams_info.append({"name": a_sam,
                          "total": total,
                          "unireads": align,
                          "percentage": "%.2f%%" % (usable_percentage * 100)})

    bowtie = JinjaTemplateCommand(
        name="bowtie summary",
        template=input["data_template"],
        param={"section_name": "bowtie",
               "output_sams": True,
               "sams": sams_info})

    write_into(bowtie, output["sum_section"])

def peaks_parse(input):
    total = 0
    fc20n = 0
    fc10n = 0
    peaks_info = {}
    with open(input) as peaks_xls:
        for line in peaks_xls:
            if line.startswith('# tags after filtering in treatment'):
                # tags after filtering in treatment: 13438948
                peaks_info["uniloc"] = int(line.strip().split(":")[1])
            if line.startswith('# d'):
                peaks_info["distance"] = int(line.strip().split("=")[1])
            if line.strip() != "" and not line.startswith("#") and not line.startswith("chr\t"):
                l = line.strip().split("\t")
                total += 1
                ## column 7th denotes fold change value
                fc = float(l[7])
                if fc >= 20:
                    fc20n += 1
                if fc >= 10:
                    fc10n += 1
    peaks_info["totalpeak"] = total
    peaks_info["peaksge20"] = fc20n
    peaks_info["peaksge10"] = fc10n
    peaks_info["peaksge20ratio"] = peaks_info["peaksge20"] / peaks_info["totalpeak"]
    peaks_info["peaksge10ratio"] = peaks_info["peaksge10"] / peaks_info["totalpeak"]
    return peaks_info

def macs2_summary(input = {"data_template": "", "macs2_peaks_xls": ""},
                  output = {"sum_section": ""}, param = {}):
    """ unique location, total peaks, peaks overlap with DHS sites, DHS sites ratio
    peaks with fold change >= 20(and ratio),
    peaks with fold change >= 10(and ratio),
    distance = shiftsize * 2
    velcro peaks number( and ratio )
    peaks = {'name1':a, 'total1': 5...} **args
    """
    peaks_info = peaks_parse(input["macs2_peaks_xls"])

    macs2_summary = JinjaTemplateCommand(
        name="peaks summary",
        template=input["data_template"],
        param={"section_name": "peaks_summary",
               "peaks": peaks_info})

    write_into(macs2_summary, output["sum_section"])

def velcro_summary(input = {"data_template": "", "macs2_peaks_xls": "", "velcro_peaks": ""}, output = {"sum_section": ""}, param={}):
    peaks_info =  peaks_parse(input["macs2_peaks_xls"])
    peaks_info["velcro"] = len(open(input["velcro_peaks"], 'r').readlines())
    peaks_info['velcropercentage'] = peaks_info["velcro"] / peaks_info["totalpeak"]

    velcro = JinjaTemplateCommand(
        name="velcro summary",
        template=input["data_template"],
        param={"section_name": "velcro",
               "peaks": peaks_info})
    write_into(velcro, output["sum_section"])

def dhs_summary(input = {"data_template": "", "macs2_peaks_xls": "", "dhs_peaks": ""}, output = {"sum_section": ""}, param={}):
    peaks_info =  peaks_parse(input["macs2_peaks_xls"])
    peaks_info["dhs"] = len(open(input["dhs_peaks"], 'r').readlines())
    peaks_info['dhspercentage'] = peaks_info["dhs"] / peaks_info["totalpeak"]
    DHS = JinjaTemplateCommand(
        name="peaks summary",
        template=input["data_template"],
        param={"section_name": "dhs",
               "peaks": peaks_info})
    write_into(DHS, output["sum_section"])

# def file_summary(input = "", output ="", param = {"has_reps": True, "files": ""}):
#     """

#     Arguments:
#     - `input`:
#     - `output`:
#     - `param`:
#     """
#     return
