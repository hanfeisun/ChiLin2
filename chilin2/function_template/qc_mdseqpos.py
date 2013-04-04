import re
import json
from chilin2.jinja_template_render import JinjaTemplateCommand, write_into

def _extract_traverse_tree(tree):
    if tree.get("node", None):
        return []

    result_list = [tree['node']]

    for a_child in tree['children']:
        result_list.extend(_extract_traverse_tree(a_child))
    return result_list

# TODO: parse ?
def qc_mdseqpos_parse_and_filter_by_z_score(input = {"latex_template": "", "seqpos": ""}, output={"latex_section": ""}, param = {"z_score_cutoff":-15}):
    """parrase mdsepose html file"""
    z_score_cutoff = param["z_score_cutoff"]
    seqpos_html_content = open(input['seqpos']).read()
    motif_tree_json_content = re.findall(
        r'var mtree = (.*)',
        seqpos_html_content)[0]

    motif_tree_json_content = motif_tree_json_content.replace("\'", "\"")
    motif_tree_dict = json.loads(motif_tree_json_content)
    all_motif_list = _extract_traverse_tree(motif_tree_dict)
    satisfied_motif_list = []
    satisfied_count = 0

    for a_motif in all_motif_list:
        if a_motif['zscore'] == 'None':
            a_motif['zscore'] = 65535
        if a_motif['factors'] == []:
            a_motif['factors'] = ['denovo']
    all_motif_list.sort(key=lambda x:x['zscore'])

    for a_motif in all_motif_list:

        if a_motif['id'].find('observed')>0:
            continue
        if satisfied_count == 10:
            break

        # z_score is a negative score, the smaller, the better
        if a_motif['zscore'] < z_score_cutoff :
            satisfied_count += 1
            satisfied_motif_list.append(a_motif)

    if satisfied_count < 5:
        satisfied_motif_list = all_motif_list[:5]
    
    motif_latex = JinjaTemplateCommand(
        name = "motif finding",
        template = input["latex_template"],
        param = {"motif_table": satisfied_motif_list,
                 "section_name": "motif"}
        )
    write_into(motif_latex, output["latex_section"])

    return satisfied_motif_list

def end_tex(input = "", output = {"latex_section": ""}, param = {}):
    end_latex = JinjaTemplateCommand(
        name = "end of latex document",
        template = input,
        param = {"section_name": "ending"})
    
    write_into(end_latex, output["latex_section"])

