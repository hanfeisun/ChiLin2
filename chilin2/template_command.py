"""
separate latex template to several object
"""
from jinja2 import Environment, FileSystemLoader
from chilin2.command import CommandInterface

class Latex(CommandInterface):
    def __init__(self, name, template=None, tool=None, param = {}, input=[], output=[]):
        CommandInterface.__init__(self, name, template=template, tool=None, param = param, input=[], output=[])
        #ChoiceLoader
        self._env = Environment(loader = FileSystemLoader("./"), ## PackageLoader('chilin', 'template'),
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
        self.dryrender = True
    
    def get_template(self):
        self.templatefile = self._env.get_template(self._template)
        
    def render(self):
        self.get_template()
        print(self.templatefile.render(**self.param))
        return True
    
    def real_invoke(self):
        """ Real-run current command"""
        print(self.param)
        self.render()
    
    def invoke(self):
        """ Invoke the command """
        self.before_invoke()
        if self._dry_run:
            self.result = self.dry_invoke()
        else:
            self.result = self.real_invoke()
        self.after_invoke()
        
    def dry_invoke(self):
        """ Dry-run current command: Pretend to run but not invoke anything """
        print("Rendering Latex part %s" % self.name)
        
    @property
    def env(self):
        return self._env
    
    @env.setter
    def env(self, e): ## setter to rule
        assert isinstance(e, Environment)
        self._env = e
    
    @property
    def err(self):
        return self.param
    
    @err.setter
    def err(self, e):
        """ if error, still render pdfs """
        self.param['end'] = e