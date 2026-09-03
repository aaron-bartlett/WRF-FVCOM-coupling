"""
Class to instantiate NU-WRF regression tests.
"""
from __future__ import print_function


class RegTest(object):
    def __init__(self, name):
        self.name = name
        self.compiler = "intel"
        self.npes = list()
        self.expected_output = ""
        self.wall_time = "2:00:00"
        self.qos = "allnccs"
        self.verification = "run"
        self.repo_type = "git"
        self.repo_url = ""
        self.repo_branch = "master"
        self.scratch_dir = "."
        self.model_dir = "."
        self.exe_dir = ""
        self.scripts_dir = "."
        self.data_dir = "."
        self.baseline_dir = "."
        self.mail_to = ""
        self.sponsor_id = ""
        self.use_batch = "no"
        self.use_modules = "no"
        self.make_system = "makeOld"
        self.use_html = "yes"
        self.clean_scratch = "yes"
        self.components = ""
        self.results = [self.name, "    -     ", "    -     ", "    -     "]
        self.debug = False
        self.update_base = "no"
        self.time_stamp = ""

    def get_opt(self, key):
        if key == "compiler":
            return self.compiler
        elif key == "npes":
            return self.npes
        elif key == "wall_time":
            return self.wall_time
        elif key == "qos":
            return self.qos
        elif key == "verification":
            return self.verification
        elif key == "repo_type":
            return self.repo_type
        elif key == "repo_url":
            return self.repo_url
        elif key == "repo_branch":
            return self.repo_branch
        elif key == "scratch_dir":
            return self.scratch_dir
        elif key == "exe_dir":
            return self.exe_dir
        elif key == "model_dir":
            return self.model_dir
        elif key == "scripts_dir":
            return self.scripts_dir
        elif key == "data_dir":
            return self.data_dir
        elif key == "baseline_dir":
            return self.baseline_dir
        elif key == "mail_to":
            return self.mail_to
        elif key == "sponsor_id":
            return self.sponsor_id
        elif key == "use_batch":
            return self.use_batch
        elif key == "use_modules":
            return self.use_modules
        elif key == "make_system":
            return self.make_system
        elif key == "use_html":
            return self.use_html
        elif key == "clean_scratch":
            return self.clean_scratch
        elif key == "time_stamp":
            return self.time_stamp
        elif key == "components":
            return self.components
        elif key == "expected_output":
            return self.expected_output
        elif key == "debug":
            return RegTest.str2bool(self.debug)
        else:
            return None
        
    @staticmethod
    def str2bool(v):
        return v.lower() in ("yes", "true", "t", "1")
    
    def get_npes_val(self, comp_index):
        return self.npes[comp_index]

    def set_opts(self, test_config, section_name):
        for key, value in test_config.items(section_name):
            if key == "compiler":
                self.compiler = value
            elif key == "npes":
                items = value.split(",")
                self.npes = [n for n in items]
            elif key == "wall_time":
                self.wall_time = value
            elif key == "qos":
                self.qos = value
            elif key == "verification":
                self.verification = value
            elif key == "repo_type":
                self.repo_type = value
            elif key == "repo_url":
                self.repo_url = value
            elif key == "repo_branch":
                self.repo_branch = value
            elif key == "scratch_dir":
                self.scratch_dir = value
            elif key == "exe_dir":
                self.exe_dir = value
            elif key == "model_dir":
                self.model_dir = value
            elif key == "scripts_dir":
                self.scripts_dir = value
            elif key == "data_dir":
                self.data_dir = value
            elif key == "baseline_dir":
                self.baseline_dir = value
            elif key == "mail_to":
                self.mail_to = value
            elif key == "sponsor_id":
                self.sponsor_id = value
            elif key == "use_batch":
                self.use_batch = value
            elif key == "use_modules":
                self.use_modules = value
            elif key == "make_system":
                self.make_system = value
            elif key == "use_html":
                self.use_html = value
            elif key == "clean_scratch":
                self.clean_scratch = value
            elif key == "time_stamp":
                self.time_stamp = value
            elif key == "components":
                self.components = value
            elif key == "expected_output":
                self.expected_output = value
            elif key == "debug":
                self.debug = RegTest.str2bool(value)
            elif key == "update_base":
                self.update_base = value

    def dump_all(self):
        for attr in dir(self):
            if hasattr(self, attr):
                print("obj.%s = %s" % (attr, getattr(self, attr)))

    def dump(self):
        options = vars(self)
        print(", ".join("%s: %s" % item for item in list(options.items())))

    def test_to_parser(self):
        try:
            import ConfigParser

            parser = ConfigParser.RawConfigParser()
        except ImportError:
            from configparser import RawConfigParser

            parser = RawConfigParser()
        parser.add_section(self.name)
        attributes = vars(self)
        for k, v in list(attributes.items()):
            parser.set(self.name, k, v)
        return parser

    def parser_to_test(self, parser):
        self.set_opts(parser, self.name)
