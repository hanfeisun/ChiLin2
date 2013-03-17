import os
import copy
from pprint import pprint
class AbstractCommand(object):
    def __init__(self, name, template=None, tool=None, param = {},input=[], output=[] ):
        self.name = name
        self.input = input
        self.output = output
        self.param = param
        self.template = template
        self.tool = tool
        self.result = None
        self.verbose_level = 1
        self.dry_run_mode = False
        self.resume = False

        # if the parent of a command is itself, it's a root
        self._parent = self
        self._commands = []

    def add(self, command):
        """ For workflow only, add a command into current workflow """
        raise NotImplemented

    def update(self, **kwargs):
        for k, v in kwargs.items():
            getattr(self, k).update(v)
        return self

    def set(self, **kwargs):
        for k,v in kwargs.items():
            setattr(self, k, v)
        return self

    def invoke(self):
        """
        Invoke the command, return True if no trouble encountered.
        Dry run mode check files `dangling` for only input
        Non-dry run mode check `missing` for both input and output
        """
        dangling_inputs = self._dangling_inputs
        if dangling_inputs:
            print("These input files might be dangling in ", self.name, dangling_inputs)
            return False

        if self.dry_run_mode:
            self.result = self._simulate()
            return True

        missing_i = self._missing_inputs
        if missing_i:
            print("These input files might be missing in", self.name, missing_i)
            return False
        self.result = self._execute()
        missing_o = self._missing_output
        if missing_o:
            print("These output files might be missing in", self.name, missing_o)
            return False
        return True



    def set_option(self, verbose_level=1, dry_run=False, recursive=True, resume = False):
        """
        :type verbose_level: int
        verbose_level:  1 - only show fatal errors          (quiet mode)
                        2 - show fatal errors and warnings  (normal mode)
                        3 - show workflow details           (verbose mode)
                        4 - debug mode                      (debug mode)
        """
        self.verbose_level = verbose_level
        self.dry_run_mode = dry_run
        self.resume = resume

    @property
    def clone(self):
        return copy.deepcopy(self)

    @property
    def _is_root(self):
        return self._parent == self

    def _simulate(self):
        """ Hook method for `invoke` in dry run mode: Pretend to run but not invoke anything """
        pass

    def _execute(self):
        """
        Hook method for `invoke`
        Return True if no trouble encountered
        """
        return True

    @property
    def _dummy_files(self):
        """
        Return a list of files that produced by leaves before current node in the tree
        """
        if self._is_root:
            # if current command is not a leaf of a tree, it shouldn't have dummy files.
            return []
        ret = []
        for a_command in self._root:
            if a_command == self:
                break
            else:
                ret += a_command._outputs
        return ret

    def _missing(self, files):
        missing = []
        for i in files:
            if not os.path.exists(i):
                missing.append(i)
                continue
            if os.path.isfile(i) and os.path.getsize(i) == 0:
                missing.append(i)
                continue
        return missing

    @property
    def _missing_inputs(self):
        """
        Return a list of files that are current command's input but doesn't exist in filesystem .
        This method is called before real invoke current command.
        """
        return self._missing(self._inputs)

    @property
    def _missing_outputs(self):
        """
        Return a list of files that are current command's output but doesn't exist in filesystem .
        This method is called after real invoke current command.
        """
        return self._missing(self._inputs)

    @property
    def _dangling_inputs(self):
        """
        Return a list of files that are current commands' input but:
        (1) doesn't exist in filesystem
        (2) can't be found as some commands' output before current command

        This Hook method is called:
        (1) on each leaf before both dry run and real dry

        If current command doesn't belong to a tree, just return missing inputs
        """
        if self._is_root:
            return self._missing_inputs
        else:
            return [i for i in self._missing_inputs if i not in self._dummy_files]

    @property
    def _root(self):
        """
        Return the root of the tree in which current command is
        """
        return self if self._is_root else self._parent._root

    def _collect(self, obj):
        if isinstance(obj, dict):
            return obj.values()
        elif isinstance(obj, str):
            return [obj]
        else:
            return obj

    @property
    def _inputs(self):
        """ Return the inputs as a list """
        return self._collect(self.input)

    @property
    def _outputs(self):
        """ Return the outputs as a list """
        return self._collect(self.output)


class ShellCommand(AbstractCommand):
    def _simulate(self):
        print("Dry invoke %s\t<=\t"%self.name, self._render())

    def _execute(self):
        print("Executing: ", self._render())

        can_skip = not self._missing_outputs
        if self.resume and can_skip:
            print("Resumed from existing result.. Skip")
            return True
        else:
            return os.system(self._render())

    def _render(self):
        """ Method that return the rendered content  """
        cmd = self.template.format(input=self.input,output=self.output,param=self.param,tool=self.tool)
        return cmd

class PythonCommand(AbstractCommand):
    def _render(self):
        return "%s < %s > %s" %(self.template, self._inputs, self._outputs)
    def _execute(self):
        print("Executing Function: ", self._render())
        return self.template(input=self.input, output=self.output, param = self.param)
    def _simulate(self):
        print("Dry invoke %s\t<=\t"%self.name, self._render() )
        return None



class Workflow(AbstractCommand):
    def __init__(self, name, *args):
        AbstractCommand.__init__(self, name, template=None, tool=None, param = {},input=[], output=[])
        self._commands = []

    def __iter__(self):
        for cmd_i in self._commands:
            if isinstance(cmd_i, Workflow):
                for cmd_i_j in cmd_i:
                    yield cmd_i_j
            else:
                yield cmd_i

    def add(self, command):
        command._parent = self
        self._commands.append(command)
        return self

    @property
    def _dangling_inputs(self):
        dangling_dict = {}
        for cmd in self:
            not_dangling = not cmd._dangling_inputs
            if not_dangling:
                continue
            if cmd.name in dangling_dict:
                # update by union of two set
                dangling_dict[cmd.name] |= set(not_dangling)
            else:
                dangling_dict[cmd.name] = set(not_dangling)
        return dangling_dict

    def invoke(self):
        dangling_inputs = self._dangling_inputs
        if dangling_inputs:
            print("The following files might be dangling. Please check whether they exists")
            pprint(dangling_inputs)
            return False
        for cmd in self:
            success_invoked = cmd.invoke()
            if not success_invoked:
                return False
        return True

    def set_option(self, verbose_level=1, dry_run=False, recursive=True, resume=False):
        self.verbose_level = verbose_level
        self.dry_run_mode = dry_run
        self.resume = resume
        if recursive:
            for cmd in self._commands:
                cmd.set_option(verbose_level, dry_run, True, resume)





