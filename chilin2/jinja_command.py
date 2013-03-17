"""
separate latex template to several object
"""
from jinja2 import Environment, FileSystemLoader
from chilin2.command import AbstractCommand

class JinjaCommand(AbstractCommand):
    def __init__(self, name, template=None, tool=None, param = {}, input=[], output=[]):
        AbstractCommand.__init__(self, name, template=template, tool=None, param = param, input=[], output=[])
        self.env = Environment(loader = FileSystemLoader("./"),
                        block_start_string = '\BLOCK{',
                        block_end_string = '}',
                        variable_start_string = '\VAR{',
                        variable_end_string = '}',
                        comment_start_string = '\#{',
                        comment_end_string = '}',
                        line_statement_prefix = '%-',
                        line_comment_prefix = '%#',
                        trim_blocks = True,
                        autoescape = False,
                        )
        self._t = self.env.get_template(self.template)

    def _execute(self):
        """ Real-run current command"""
        self.result = self._t._render(**self.param)

        
    def _simulate(self):
        """ Dry-run current command: Pretend to run but not invoke anything """
        print("Rendering Latex part %s" % self.name)
