"""
separate latex template to several object
"""
import subprocess
from jinja2 import Environment, FileSystemLoader, PackageLoader
from pyflow.command import AbstractCommand

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

def write_into(jinja_template, file_path):
    jinja_template.invoke()
    with open(file_path, "w") as f:
        f.write(jinja_template.result)

def write_and_run_Rscript(jinja_template_r, file_path):
    write_into(jinja_template_r, file_path)
    subprocess.call("Rscript %s" %file_path, shell=True)



