from collections import OrderedDict
from os import environ, path, getenv
import json
import logging

logger = logging.getLogger("NuwrfBuild")

class NuwrfBuild(object):
    version = "v11 Ekman"
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
        "wrfonly",
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
        "envs", # used to "dump" environment
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

    all_targets = dict(target_dir)
    del all_targets["met"]

    allchem_targets = dict(all_targets)
    del allchem_targets["wrf"]
    allchem_targets["chem"] = "WRF"
    
    allkpp_targets = dict(all_targets)
    del allkpp_targets["wrf"]
    allkpp_targets["kpp"] = "WRF"
    
    def __init__(self):
        self.build_config = None
        self.targets = list()
        self.target_options = dict()
        self.options = list()
        self.dump_envs = False
        self.prefix = None
        self.config = None
        self.ideal_case = False
        self.old_build_state = dict()
        self.build_settings_file = NuwrfBuild.top_dir + "/.build_settings"
        self.build_state = { "version":NuwrfBuild.version,
                             "libdir_tag":None,
                             "prefix":self.prefix,
                             "config_file":self.config,
                             "debug_opt":"0",
                             "nompi_opt": "0",
                             "ideal_case_opt":"0",
                             "chem_opt":"0",
                             "kpp_opt":"0",
                             "wrflis_opt":"0",
                             "nests_opt":"1",
                             "wrf_configure_opt":None,
                             "wps_configure_opt":None,
                             "wrf_configure_lis":None,
                             "components_built":list(),
                             "modules_used":list()
                             }

    @staticmethod
    def set_env_vars():
        for k, v in environ.items():
            NuwrfBuild.env_vars[k] = v

    @staticmethod
    def dump_env_vars():
        with open(NuwrfBuild.top_dir + "/nu-wrf.envs", "w") as f:
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
            if "WRF_CHEM" in NuwrfBuild.env_vars or "lis" in name:
                self.target_options[name].remove("nest=2")
                self.target_options[name].remove("nest=3")

    def create_build_config(self):
        logger.info("Set .build_settings file with NU-WRF build options.")
        self.build_state["LIBDIR_TAG"] = NuwrfBuild.env_vars["LIBDIR_TAG"]
        self.build_state["config_file"] = self.config
        self.build_state["install_dir"] = self.prefix
        self.build_state["wrf_configure_opt"] = NuwrfBuild.env_vars["WRF_CONFIGURE_MPI_OPT"]
        self.build_state["wps_configure_opt"] = NuwrfBuild.env_vars["WPS_CONFIGURE_MPI_OPT"]
        self.build_state["wrf_configure_lis"] = NuwrfBuild.env_vars["WRF_CONFIGURE_LIS_MPI"]
        try:
            with open(NuwrfBuild.top_dir + "/.modules") as f:
                mods_used = f.readlines()[0].rstrip().split(" ")
                self.build_state["modules_used"] = mods_used[4::2]
        except IOError:
            self.build_state["modules_used"] = " NONE"
        self.set_wrf_envs()
        self.update_build_config_opts()
        self.write_build_config()


    def read_build_config(self):
        def is_a_in_b_(a, b):
            """ return True if list a is in b """
            return set(set(a)).issubset(set(b))

        with open(self.build_settings_file) as f:
            self.old_build_state = json.load(f)
        self.build_state["version"] = NuwrfBuild.version
        self.build_state["libdir_tag"] = getenv("LIBDIR_TAG")
        self.build_state["prefix"] = self.old_build_state["prefix"]
        self.build_state["config_file"] = self.old_build_state["config_file"]
        self.build_state["debug_opt"] = self.old_build_state["debug_opt"]
        self.build_state["nompi_opt"] = self.old_build_state["nompi_opt"]
        self.build_state["nests_opt"] = self.old_build_state["nests_opt"]
        self.build_state["chem_opt"] = self.old_build_state["chem_opt"]
        self.build_state["kpp_opt"] = self.old_build_state["kpp_opt"]
        self.build_state["ideal_case_opt"] = self.old_build_state["ideal_case_opt"]
        self.build_state["wrflis_opt"] = self.old_build_state["wrflis_opt"]
        self.build_state["wrf_configure_opt"] = self.old_build_state["wrf_configure_opt"]
        self.build_state["wps_configure_opt"] = self.old_build_state["wps_configure_opt"]
        self.build_state["wrf_configure_lis"] = self.old_build_state["wrf_configure_lis"]

        try:
            with open(NuwrfBuild.top_dir + "/.modules") as f:
                mods_used = f.readlines()[0].rstrip().split(" ")
                self.build_state["modules_used"] = mods_used[4::2]
        except IOError:
            self.build_state["modules_used"] = " NONE"
        self.set_wrf_envs()
        self.update_build_config()

        self.build_state["components_built"] = self.old_build_state.get("components_built")
        if not is_a_in_b_(self.build_state["modules_used"],
                          self.old_build_state.get("modules_used")):
            logger.error("Module environment has changed. Run './build.sh allclean' first.")
            print("Curr:",self.build_state["modules_used"])
            print("Prev:",self.old_build_state.get("modules_used"))
            import sys
            sys.exit()
        else:
            self.build_state["modules_used"] = self.old_build_state.get("modules_used")

    def update_built_components(self):
        built_values = self.build_state.get("components_built")
        for t in self.targets:
            built_values.append(t)
        foo = list(set(built_values))
        self.build_state["components_built"] = foo

    def update_build_config_opts(self):
        if "debug" in self.options:
            self.build_state["debug_opt"] = "1"

        if "nompi" in self.options:
            self.build_state["nompi_opt"] = "1"

        if "chem" in self.targets:
            self.build_state["chem_opt"] = "1"

        if "kpp" in self.options:
            self.build_state["kpp_opt"] = "1"

        if "wrflis" in self.options:
            self.build_state["wrflis_opt"] = "1"

        if "ideal_case" in self.options:
            self.build_state["ideal_case_opt"] = "1"

        nestx = [s for s in self.options if "nest" in s]
        if len(nestx) == 1:
            if "1" in nestx and self.old_build_state["nests_opt"] != "1":
                nests = nestx[0].split("=")[-1]
                self.build_state["nests_opt"] = nests

    def update_build_config(self):
        clean_first = "cleanfirst" in self.options

        if "debug" in self.options and self.old_build_state["debug_opt"] == "0":
            self.build_state["debug_opt"] = "1"
        else:
            self.build_state["debug_opt"] = "0"

        if "nompi" in self.options and self.old_build_state["nompi_opt"] == "0":
            self.build_state["nompi_opt"] = "1"
        else:
            self.build_state["debug_opt"] = "0"

        if "chem" in self.options and self.old_build_state["chem_opt"] == "0":
            self.build_state["chem_opt"] = "1"
        else:
            self.build_state["chem_opt"] = "0"

        if "kpp" in self.options and self.old_build_state["kpp_opt"] == "0":
            self.build_state["kpp_opt"] = "1"
        else:
            self.build_state["kpp_opt"] = "0"

        if "wrflis" in self.options and self.old_build_state["wrflis_opt"] == "0":
            self.build_state["wrflis_opt"] = "1"
        else:
            self.build_state["wrflis_opt"] = "0"

        if "ideal_case" in self.options and self.old_build_state["ideal_case_opt"] == "0":
            self.build_state["ideal_case_opt"] = "1"
        else:
            self.build_state["ideal_case_opt"] = "0"

        nestx = [s for s in self.options if "nest" in s]
        if len(nestx) == 1:
            if "1" in nestx and self.old_build_state["nests_opt"] != "1":
                nests = nestx[0].split("=")[-1]
                self.build_state["nests_opt"] = nests
            else:
                self.build_state["nests_opt"] = "1"


    def remove_built_component_opt(self, item):
        try:
            self.build_state["components_built"].remove(item)
        except ValueError:
            logger.debug("Will rebuild "+item)

    def write_build_config(self):
        foo = json.dumps(self.build_state)
        with open(self.build_settings_file, "w") as f:
            f.write(foo)

    def dump_build_settings(self):
        print("Current build settings")
        print("----------------------")
        for key, val in self.build_state.items():
            print("{0} = {1}".format(key, val))

    def set_wrf_envs(self):
        envs = NuwrfBuild.env_vars
        # WRF-LIS coupling is the default build option
        envs["WRF_LIS"] = "1"
        if "wrfonly" in self.targets:
            envs["WRF_LIS"] = "0"
        if "WRF_CHEM" not in envs:
            envs["WRF_CHEM"] = "0"
        if "WRF_KPP" not in envs:
            envs["WRF_KPP"] = "0"
        envs["WRF_EM_CORE"] = "1"
        envs["WRF_NMM_CORE"] = "0"
        envs["WRF_NMM_NEST"] = "0"
        envs["HWRF"] = "0"
        envs["WRF_DA_CORE"] = "0"
        envs["WRF_COAMPS_CORE"] = "0"
        envs["WRF_EXP_CORE"] = "0"
        envs["WRF_TITAN"] = "0"
        envs["WRF_MARS"] = "0"
        envs["WRF_VENUS"] = "0"
        envs["OMP_NUM_THREADS"] = "1"
        envs["WRF_DFI_RADAR"] = "0"
        envs["WRF_CONVERT"] = "0"
        envs["MADIS"] = "0"
        envs["BUFR"] = "0"
        envs["RTTOV"] = "0"
        envs["CRTM"] = "0"
        envs["WRF_ESMF_LIBS"] = envs["WRF_ESMF_LIBS_MPI"]
        envs["WRF_CONFIGURE_OPT"] = envs["WRF_CONFIGURE_MPI_OPT"]
        envs["WPS_CONFIGURE_OPT"] = envs["WPS_CONFIGURE_MPI_OPT"]
        envs["WRF_CONFIGURE_LIS"] = envs["WRF_CONFIGURE_LIS_MPI"]
        if "nompi" in self.options:
            envs["WRF_ESMF_LIBS"] = envs["WRF_ESMF_LIBS_NOMPI"]
            envs["WRF_CONFIGURE_OPT"] = envs["WRF_CONFIGURE_NOMPI_OPT"]
            envs["WPS_CONFIGURE_OPT"] = envs["WPS_CONFIGURE_NOMPI_OPT"]
            envs["WRF_CONFIGURE_LIS"] = envs["WRF_CONFIGURE_LIS_NOMPI"]
