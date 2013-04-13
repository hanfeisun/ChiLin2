import re
import json
from chilin2.helpers import JinjaTemplateCommand, template_dump, json_load

def _extract_traverse_tree(tree):


    if not tree.get("node", None):
        return []

    result_list = [tree['node']]

    for a_child in tree['children']:
        result_list.extend(_extract_traverse_tree(a_child))

    return result_list


def stat_seqpos(input = {"template": "", "seqpos": ""}, output={"latex_section": ""}, param = {"z_score_cutoff":-15}):
    """parrase mdsepose html file"""
    z_score_cutoff = param["z_score_cutoff"]
    seqpos_html_content = open(input['seqpos']).read()
    js_string = re.findall(
        r'var mtree = (.*)',
        seqpos_html_content)[0]

    # TODO(hanfei) : check the issue that mdseqpos table is empty
    js_string = js_string.replace("\'", "\"")
    mdseqpos_result = json.loads(js_string)

    all_motif_list = _extract_traverse_tree(mdseqpos_result)
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

    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    result_dict["stat"]["satisfied_motifs"] = satisfied_motif_list
    with open(output["json"], "w") as f:
        json.dump(result_dict, f, indent=4)

def latex_seqpos(input, output, param):
    json_dict = json_load(input["json"])

    latex = JinjaTemplateCommand(
        name = "motif finding",
        template = input["template"],
        param = {"motif_table": json_dict["stat"]["satisfied_motifs"],
                 "section_name": "motif",
                 "render_dump": output["latex"]})

    template_dump(latex)





