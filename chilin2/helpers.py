"""
separate latex template to several object
"""
import json
from jinja2 import Environment, FileSystemLoader
from samflow.command import AbstractCommand, ShellCommand

env = Environment(loader = FileSystemLoader("/"),
    block_start_string = '\BLOCK{',
    block_end_string = '}',
    variable_start_string = '\VAR{',
    variable_end_string = '}',
    comment_start_string = '\#{',
    comment_end_string = '}',
    line_statement_prefix = '%-',
    line_comment_prefix = '%#',
    trim_blocks = True,
    autoescape = False,)

def surround_by_quote(a_list):
    return ['"%s"' % an_element for an_element in a_list]

env.filters["surround_by_quote"] = surround_by_quote


class JinjaTemplateCommand(AbstractCommand):
    def __init__(self, template, tool=None, param = {}, input=[], output=[], name = ""):
        AbstractCommand.__init__(self, template=template, tool=None, param = param, input=[], output=[], name = name)
        self.env = env

        self._t = self.env.get_template(self.template)

    def _execute(self):
        """ Real-run current command"""
        self.result = self._t.render(input = self.input, output = self.output, **self.param)
        return True

    def _simulate(self):
        """ Dry-run current command: Pretend to run but not invoke anything """
        print("Rendering Latex part %s" % self.name, self.template)
        return True

def template_dump(jinja_template):
    jinja_template.invoke()
    with open(jinja_template.param["render_dump"], "w") as f:
        f.write(jinja_template.result)

def r_exec(jinja_template_r):
    ShellCommand(template="Rscript {input}",
        input=jinja_template_r.param["render_dump"],
        output=jinja_template_r.param["pdf"]).invoke()

def json_dump(json_dict):
    json_file = json_dict["output"]["json"]
    with open(json_file, "w") as f:
        json.dump(json_dict, f, indent=4)
    return json_file

def json_load(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    return json_dict

def latex_start(input = {"template": ""}, output = {"latex": ""}, param = {"id": ""}):
    end_latex = JinjaTemplateCommand(
        name = "end of latex document",
        template = input["template"],
        param = {"section_name": "begin",
                 "prefix_dataset_id": param["id"],
                 "render_dump": output["latex"]
        })

    template_dump(end_latex)
def latex_end(input = {"template": ""}, output = {"latex": ""}, param = {}):
    end_latex = JinjaTemplateCommand(
        name = "end of latex document",
        template = input["template"],
        param = {"section_name": "ending",
                 "render_dump": output["latex"]
                 })

    template_dump(end_latex)


def underline_to_space(x):
    if type(x) == str:
        return x.replace("_", " ")
    return x