import jinja2
from chilin2.jinja_template_render import JinjaTemplateCommand, write_into

def bowtie_summary(input = {"data_template": ""}, output = {"sum_section": ""},
                   param = {"bowtie_summary": "", "peaks_xls": "", "sam_files": ""}):
    """bowtie data summary
    sams = [{'name1':a, 'total1': 5...}, {'name2':c, 'total2': 3...}...] **args
    """
    ## TODO: unique location ?
    sams = []
    for i in range(len(param["bowtie_summary"])):
        ## bowtie summary has same length with macs2 xls
        # with open(param["peaks_xls"][i]) as macs:
        #     content = macs.readlines()
        #     uniloc = int(content[16].strip().split(":")[1])
        with open(param["bowtie_summary"][i]) as bowtie_sum:
            content = bowtie_sum.readlines()
            total = int(content[0].strip().split(":")[1])
            align = int(content[4].strip().split()[1])
            usable_percentage = float(align)/total
        sams.append({"name": input["sam_files"][i],  "total": total,
                     "unireads": align, "percentage": str(usable_percentage*100) + "%"})

    print(sams)
    bowtie = JinjaTemplateCommand(
        name = "bowtie summary",
        template = input["data_template"],
        param = {"section_name": "bowtie",
                 "output_sams": True,
                 "sams": sams})
    write_into(bowtie ,output["sum_section"])

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
