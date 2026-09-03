from collections import OrderedDict
from os import environ, path
from sys import version_info

try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser


class NuwrfBuild(object):
    
    version = "v10"
    env_vars = OrderedDict()
    top_dir = environ.get("NUWRFDIR")
    
    # A list of all valid NU-WRF build targets:
    valid_targets = [
        "all",
        "allclean",
        "allchem",
        "allkpp",
        "arw",
        "chem",
        "doc",
        "gmp",
        "grad",
        "gsdsu",
        "kpp",
        "doc",
        "ldt",
        "lvt",
        "lis",
        "met",
        "rip",
        "upp",
        "utils",
        "wps",
        "wrf",
        "ideal_b_wave",
        "ideal_convrad",
        "ideal_heldsuarez",
        "ideal_les",
        "ideal_quarter_ss",
        "ideal_scm_xy",
        "ideal_scm_lis_xy",
        "ideal_tropical_cyclone",
        "casa2wrf",
        "gocart2wrf",
        "geos2wrf",
        "sst2wrf",
        "lisWrfDomain",
        "lis4scm",
        "ndviBareness4Wrf",
        "plot_chem",
        "prep_chem_sources",
    ]
    # Dictionary of top level targets and their corresponding directories
    target_dir = {
        "arw": "ARWpost",
        "doc": "docs",
        "gmp": "postproc/GMP/SRC",
        "grad": "postproc/GRAD/SRC",
        "gsdsu": "GSDSU/SRC",
        "lis": "LISF/lis",
        "ldt": "LISF/ldt",
        "lvt": "LISF/lvt",
        "met": "MET",
        "rip": "RIP4",
        "upp": "UPP",
        "utils": "utils",
        "wps": "WPS",
        "wrf": "WRF",
    }
    utils_exe = [
        "casa2wrf",
        "geos2wrf",
        "gocart2wrf",
        "sst2wrf",
        "ndviBareness4Wrf",
        "lisWrfDomain",
        "lis4scm",
        "plot_chem",
        "prep_chem_sources",
    ]

    def __init__(self):
        self.prefix = None
        self.config = None
        self.build_settings_file = None
        self.build_state = dict()
        self.build_config = None
        self.build_boot = True
        self.targets = list()
        self.target_options = dict()
        self.options = list()
        self.dump_envs = False

    @staticmethod
    def set_env_vars():
        for k, v in environ.items():
            NuwrfBuild.env_vars[k] = v

    @staticmethod
    def dump_env_vars():
        with open(NuwrfBuild.top_dir+"/nu-wrf.envs", "w") as f:
            f.write("# NU-WRF build configuration variables\n")
            for k, v in NuwrfBuild.env_vars.items():
                if " " in v:
                    f.write("export " + k + "=" + "'{}'".format(v) + "\n")
                else:
                    f.write("export " + k + "=" + v + "\n")

    def set_targets_options(self, target_list):
        # Set list of valid target build options
        for name in target_list:
            self.target_options.setdefault(name, [])
            self.target_options[name].append("rebuild")
            self.target_options[name].append("cleanfirst")
            self.target_options[name].append("debug")
            self.target_options[name].append("nompi")
            self.target_options[name].append("ideal_case")
            self.target_options[name].append("skip_clm4")
            self.target_options[name].append("nest=1")
            self.target_options[name].append("nest=2")
            self.target_options[name].append("nest=3")
            # WRF cannot be run coupled to LIS if preset-moves (nest=2) or
            # vortex-tracking nesting (nest=3) is used. Same for WRF-chem.
            if "BUILD_WRF_CHEM" in NuwrfBuild.env_vars or "lis" in name:
                self.target_options[name].remove("nest=2")
                self.target_options[name].remove("nest=3")

    def get_build_config(self):
        if not path.isfile(NuwrfBuild.top_dir+"/.build_settings"):
            self.build_config = self.set_build_config()
        else:
            self.build_boot = False
            self.build_config = RawConfigParser()
            self.build_config.optionxform = str
            self.build_config.read(NuwrfBuild.top_dir + "/.build_settings")
        # Loop over entries, populate build state
        for k,v in self.build_config.items("build_settings"):
            self.build_state[k] = v
            
    def set_build_config(self):
        print("Set .build_settings file with NU-WRF build options.")

        self.build_config = RawConfigParser(allow_no_value=True)
        self.build_config.optionxform = str
        self.build_config.add_section("build_settings")
        
        self.build_config.set("build_settings",
                              "NU_WRF_version: ",NuwrfBuild.version)
        self.build_config.set("build_settings",
                              "LIBDIR_TAG: ",NuwrfBuild.env_vars["LIBDIR_TAG"])
        self.build_config.set("build_settings",
                              "config_file: ",self.config)
        
        self.build_config.set("build_settings","; "+
                              "Build options:")

        self.build_config.set("build_settings", "build_debug", "0")
        if "debug" in self.options:
            self.build_config.set("build_settings", "build_debug", "1")

        self.build_config.set("build_settings", "build_ideal_case", "0")
        if "ideal_case" in self.options:
            self.build_config.set("build_settings", "build_ideal_case", "1")

        self.build_config.set("build_settings", "build_nompi", "0")
        if "nompi" in self.options:
            self.build_config.set("build_settings", "build_nompi", "1")

        self.build_config.set("build_settings", "build_chem", "0")
        if "BUILD_WRF_CHEM" in NuwrfBuild.env_vars:
            self.build_config.set("build_settings", "build_chem", "1")

        self.build_config.set("build_settings", "build_kpp", "0")
        if "BUILD_WRF_KPP" in NuwrfBuild.env_vars:
            self.build_config.set("build_settings", "build_kpp", "1")

        self.build_config.set("build_settings", "build_wrf_lis", "0")
        if "BUILD_WRF_LIS" in NuwrfBuild.env_vars:
            self.build_config.set("build_settings", "build_wrf_lis", "1")
            
        nestx = [s for s in self.options if "nest" in s]
        if len(nestx) == 1:
            nests = nestx[0].split("=")[-1]
        else:
            nests = "1"
        self.build_config.set("build_settings", "nests", nests)
        self.build_config.set("build_settings",
                                     "install_path: ",self.prefix)

        self.build_config.set("build_settings","; "+
                                     "Config file options:")

        if "WRF_CONFIGURE_MPI_OPT" in NuwrfBuild.env_vars:
            self.build_config.set(
                "build_settings",
                "wrf_configure_opt",
                NuwrfBuild.env_vars["WRF_CONFIGURE_MPI_OPT"],
            )

        if "WRF_CONFIGURE_LIS_MPI" in NuwrfBuild.env_vars:
            self.build_config.set(
                "build_settings",
                "wrf_configure_lis_file",
                NuwrfBuild.env_vars["WRF_CONFIGURE_LIS_MPI"],
            )
            
        if "WPS_CONFIGURE_MPI_OPT" in NuwrfBuild.env_vars:
            self.build_config.set(
                "build_settings",
                "wps_configure_opt",
                NuwrfBuild.env_vars["WPS_CONFIGURE_MPI_OPT"],
            )
        
        self.build_config.set("build_settings","; "+
                                     "Components built:")
        for t in self.targets:
            self.build_config.set("build_settings","build_"+t, "1")
                                         
        self.build_config.set("build_settings","; "+
                                     "Modules used:")
        try:
            with open(NuwrfBuild.top_dir+"/.modules") as f:
                mods = f.readlines()[1:]
                for line in mods:
                    self.build_config.set("build_settings","; ---  "+line.strip())
        except IOError:
            self.build_config.set("build_settings","; ---  NONE")
            pass  

        self.update_build_config()

        return self.build_config

    def get_build_config_opt(self, key):
        if key in self.build_state:
            return self.build_state[key]
        return None

    def set_build_config_opt(self, key, val):
        self.build_config[key] = val

    def remove_build_config_opt(self, key):
        if key in self.build_state:
            del self.build_state[key]

    def update_build_config(self):
        with open(NuwrfBuild.top_dir+"/.build_settings", "w") as f:
            self.build_config.write(f)

