from configparser import ConfigParser, NoSectionError
import os

class ChiLinConfig(object):
    def __init__(self, conf):
        self._verbose_level = 1
        self._conf = ConfigParser()
        self._conf.read(conf)
        self.root_dir = os.path.dirname(conf)

    def set_option(self, verbose_level=1):
        """
        :type verbose_level: int
        verbose_level:  1 - only show fatal errors          (quiet mode)
                        2 - show fatal errors and warnings  (normal mode)
                        3 - show workflow details           (verbose mode)
                        4 - debug mode                      (debug mode)
        """
        self._verbose_level = verbose_level

    def get(self, section, option):
        return self._conf.get(section, option)

    def get_path(self, section, option):
        return self.to_abs_path(self.get(section, option))

    def items(self, section):
        try:
            return self._conf.items(section)
        except NoSectionError:
            if self._verbose_level >= 2:
                print("Warning: No such section: ", section)
                print("This will return a empty dict")
            return {}

    @property
    def id(self):
        return self.get("Basis", "id")

    @property
    def target_dir(self):
        return self.get("Basis", "output")

    @property
    def target_prefix(self):
        return os.path.join(self.target_dir, self.id)

    @property
    def treatment_map(self):
        return list(zip(self._treatment_input, self._treatment_target_prefix))

    @property
    def control_map(self):
        return list(zip(self._control_input, self._control_target_prefix))

    @property
    def sample_map(self):
        return self.treatment_map + self.control_map

    def to_abs_path(self, path):
        abs_path = path
        if not os.path.isabs(path):
            abs_path = os.path.join(self.root_dir, abs_path)
        return abs_path

    @property
    def _control_input(self):
        return [self.to_abs_path(i.strip()) for i in self.get("Basis", "control").split(",")]

    @property
    def _treatment_input(self):
        return [self.to_abs_path(i.strip()) for i in self.get("Basis", "treat").split(",")]

    @property
    def _treatment_target_prefix(self):
        return [os.path.join(self.target_dir,
            self.id + "_treat_rep" + str(num+1)
        ) for num in range(len(self._treatment_input))]

    @property
    def _control_target_prefix(self):
        return [os.path.join(self.target_dir,self.id + "_control_rep" + str(num+1))
                for num in range(len(self._control_input))]



