from chilin2.jinja_template_render import JinjaTemplateCommand, write_into

def text_bowtie_summary(input={"data_template": "", "bowtie_summary": []},
                        output={"sum_section": ""},
                        param={"sam_files": []}
):
    """ sams = [{'name1':a, 'total1': 5...}, {'name2':c, 'total2': 3...}...] **args """

    # unique location is in text_macs2_summary part
    sams_info = []
    for a_summary, a_sam in zip(input["bowtie_summary"], param["sam_files"]):
        summary_lines = open(a_summary).readlines()
        print("ppyy" + a_summary)
        print(summary_lines)
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

# def macs2_summary(input = "", output = "", param = ""):
#     """
#     Arguments:
#     - `input`:
#     - `output`:
#     - `param`:
#     """
#     fhd = open(self.rule.macs.peaksxls,"r")
#     total = 0
#     fc20n = 0
#     fc10n = 0
#     for i in fhd:
#         i = i.strip()
#         if i and not i.startswith("#") and not i.startswith("chr\t"):
#             total += 1
#             fs = i.split("\t")
#             fc = fs[7]
#             if fc >= 20:
#                 fc20n += 1
#             if fc >= 10:
#                 fc10n += 1
#     self.macsinfo['totalpeak'] = total
#     self.macsinfo['peaksge20'] = fc20n
#     self.macsinfo['peaksge10'] = fc10n
#     if total != 0:
#         self.macsinfo['peaksge10ratio'] = float(fc10n)/total
#         self.macsinfo['peaksge20ratio'] = float(fc20n)/total
#     else:
#         self.macsinfo['peaksge10ratio'] = 0.0001
#         self.macsinfo['peaksge20ratio'] = 0.0001
#     self.macsinfo['distance'] = float(self.shiftsize) * 2
#     self.rendercontent['ratios'] = self.macsinfo
#     ## dhs, velcor and all peaks summary
#     lenall = len(open(self.rule.macs.treatpeaks, 'rU').readlines())
#     if a_type == 'dhs':
#         lena_type = len(open(self.rule.bedtools.dhs, 'r').readlines())
#     elif a_type == 'velcro':
#         if species:
#             lena_type = len(open(self.rule.bedtools.velcro, 'r').readlines())
#     self.ratio[a_type] = lena_type
#     if lenall != 0:
#         self.ratio[a_type + 'percentage'] = float(lena_type)/lenall
#     else:
#         self.ratio[a_type + 'percentage'] = 0.0001

#     for k, v in self.ratio.iteritems():
#         self.rendercontent['ratios'][k] = v
#     return

# def file_summary(input = "", output ="", param = {"has_reps": True, "files": ""}):
#     """

#     Arguments:
#     - `input`:
#     - `output`:
#     - `param`:
#     """
#     return
