from shutil import copyfile, move
from filecmp import cmp
import logging
import os
import shared.utils as u

logger = logging.getLogger(__name__)

wrf_dir = os.environ.get("NUWRFDIR") + "/WRF/"


def configure_lis(my_build):
    envs = my_build.env_vars
    options = my_build.options
    if os.environ.get("DEBUG_BUILD") is not None:
        return 0
    opt_para = ""
    if "nompi" in options:
        opt_para = "0"
        envs["WRF_ESMF_LIBS"] = envs["WRF_ESMF_LIBS_NOMPI"]
    opt_optim = ""
    if "debug" in options:
        opt_optim = "-1"
        
    os.chdir(os.environ.get("NUWRFDIR") + "/" + my_build.target_dir["lis"])

    heredoc = opt_para+"\n"+opt_optim+"\n"+"\n"*15+"\nEOF"
    rc = u.run_configure(r"./configure <<EOF ", heredoc, envs)
    logger.debug("rc=" + str(rc))
    
    os.chdir("make")
    if os.environ.get("DEBUG_BUILD") is not None:
        return 0
    if not os.path.isfile("configure.lis"):
        return 1

    # NOMPI:
    with open ("configure.lis") as f:
        lines = f.readlines()
    if "nompi" in options:
        lnum = 0
        with open("configure.lis","w") as f:
            for line in lines:
                lnum += 1
                if lnum < 4:
                    continue
                elif lnum == 4:
                    f.write("FC = " + envs["LIS_FC_NOMPI"] + "\n")
                    f.write("FC77 = " + envs["LIS_FC_NOMPI"] + "\n")
                    f.write("LD = " + envs["LIS_FC_NOMPI"] + "\n")
                    f.write("CC = " + envs["LIS_CC_NOMPI"] + "\n")
                elif lnum > 4:
                    if line.startswith("FFLAGS"):
                        line = line.strip()+" -DHIDE_MPI\n"
                    f.write(line)

        u.sed_inplace(
            "configure.lis",
            r"^MOD_ESMF.*",
            "MOD_ESMF = " + envs["ESMF_NOMPI_MOD"],
        )

        u.sed_inplace(
            "configure.lis",
            r"^LIB_ESMF.*",
            "LIB_ESMF = " + envs["ESMF_NOMPI_LIB"],
        )
            
    # Update Filepath and LIS_plugins
    logger.debug("Update Filepath and LIS_plugins")
    u.sp_call("python plugins.py")
    os.chdir(os.environ.get("NUWRFDIR") + "/" + my_build.target_dir["lis"])
    return rc              
    
def configure_wrf(my_build):
    envs = my_build.env_vars
    options = my_build.options
    # Set nest option for WRF. Default is nest=1:
    nest = "1"
    if "nest=2" in options:
        nest = "2"
    elif "nest=3" in options:
        nest = "3"
    # Reset nest option for NOMPI case (ideal cases)
    if "nompi" in options:
        nest = "0"
        
    os.chdir(os.environ.get("NUWRFDIR") + "/" + my_build.target_dir["wrf"])
    logger.debug("Generate configure.wrf")
    heredoc = "<<EOF\n" + envs["WRF_CONFIGURE_OPT"] + "\n" + nest + "\nEOF"
    if "debug" in options:
        rc = u.run_configure("./configure -D ", heredoc, envs)
    else:
        rc = u.run_configure("./configure -s ", heredoc, envs)
    logger.debug("rc=" + str(rc))

    if os.environ.get("DEBUG_BUILD") is not None:
        return 0
    if not os.path.isfile("configure.wrf"):
        return 1

    logger.debug("Edit configure.wrf")
    # Optionally disable compilation of CLM4 in WRF
    if "WRF_SKIP_CLM4" in envs:
        u.sed_inplace("configure.wrf", r"-DWRF_USE_CLM", "")

    # Remove NETCDF4 preprocessor flag if appropriate
    if "NETCDF4" not in envs:
        u.sed_inplace("configure.wrf", r"-DUSE_NETCDF4_FEATURES", "")
        
    # FVCOM coupling
    fvcom_cpl_opt = "-DFVCOM_CPL"  # Default CPP flag
    if envs["FVCOM_CPL"] == "0":
        u.replace_infile("configure.wrf", fvcom_cpl_opt, "")

    fvcom_cpl_opt = "FVCOM_CPL\t=\t0" # Option used by Config.pl
    if envs["FVCOM_CPL"] == "1":
        u.replace_infile("configure.wrf", fvcom_cpl_opt, "FVCOM_CPL\t=\t1")

    # WRF electrification scheme hacks
    wrf_elec_opt = "-DWRF_ELEC"  # Default CPP flag
    if envs["WRF_ELEC"] == "0":
        u.replace_infile("configure.wrf", wrf_elec_opt, "")

    wrf_elec_opt = "WRF_ELEC\t=\t0" # Option used by Config.pl
    if envs["WRF_ELEC"] == "1":
        u.replace_infile("configure.wrf", wrf_elec_opt, "WRF_ELEC\t=\t1")

    # Build with WRF-LIS coupling when wrf,chem targets are selected
    # Otherwise, clear WRF_LIS flags when wrfonly target is selected
    wrf_lis_opt = "-DWRF_LIS"  # Default CPP flag
    if envs["WRF_LIS"] == "0":
        u.replace_infile("configure.wrf", wrf_lis_opt, "")

    wrf_lis_opt = "WRF_LIS\t=\t1" # Option used by Config.pl
    if envs["WRF_LIS"] == "0":
        u.replace_infile("configure.wrf", wrf_lis_opt, "WRF_LIS\t=\t0")

    if "debug" in options:
        u.sed_inplace(
            "configure.wrf",
            r"^CFLAGS_LOCAL.*",
            "CFLAGS_LOCAL = " + envs["WRF_DEBUG_CFLAGS_LOCAL"],
        )
        u.sed_inplace(
            "configure.wrf", r"^FCOPTIM.*", "FCOPTIM = " + envs["WRF_DEBUG_FCOPTIM"]
        )
        u.sed_inplace(
            "configure.wrf", r"^FCNOOPT.*", "FCNOOPT = " + envs["WRF_DEBUG_FCNOOPT"]
        )
        u.sed_inplace(
            "configure.wrf",
            r"^FCREDUCEOPT.*",
            "FCREDUCEOPT = " + envs["WRF_DEBUG_FCNOOPT"],
        )
    return rc


def configure(my_build, target):
    logger.info("Configure [" + target + "]")
    envs = my_build.env_vars
    options = my_build.options

    if target == "arw":
        heredoc = "<<EOF\n " + envs["ARWPOST_CONFIGURE_OPT"] + " \nEOF"
        rc = u.run_configure("./configure ", heredoc, envs)

    elif target == "lvt":
        opt_level = ""
        if any("debug" in x for x in options):
            opt_level = "-1"
        # A very ugly and temporary workaround for LVT
        envs["LD_LIBRARY_PATH"] = envs["LVT_HDF5"]+"/lib:"+os.getenv("LD_LIBRARY_PATH") 
        heredoc = opt_level+"\n" + "\n"*13+"\nEOF"
        rc = u.run_configure(r'./configure <<EOF ', heredoc, envs)

    elif target == "ldt":
        o_lvl = "2"
        if any("debug" in x for x in options):
            o_lvl = "-1"
        o_geotif = "0"
        o_lgeotif = "0"
        if envs["LDT_USE_GEOTIFF"] == "1":
            o_geotif = "1"
            o_lgeotif = "1"
        o_hist = "1"
        if envs["LDT_HIST_TIMESTAMP"] == "0":
            o_hist = "0"
        
        heredoc = "\n"+o_lvl+"\n"+"\n"*9+o_geotif+"\n"+o_lgeotif+"\n"+o_hist+"\n"+"\nEOF"
        rc = u.run_configure(r'./configure <<EOF ', heredoc, envs)

    elif target == "met":
        u.sp_call("touch config.h.in")
        u.sp_call("touch aclocal.m4")
        fs = u.find_files(".", "Makefile.am")
        for f in fs:
            u.sp_call("touch " + f)
        fs = u.find_files(".", "Makefile.in")
        for f in fs:
            u.sp_call("touch " + f)
        opts = (
            "--prefix="
            + os.getcwd()
            + " --enable-grib2 --enable-modis --enable-mode_graphics --enable-lidar2nc "
        )
        rc = u.run_shell_command("./configure " + opts, envs)

    elif target == "rip":
        heredoc = "<<EOF \n" + envs["RIP_CONFIGURE_OPT"] + "\nEOF"
        rc = u.run_configure("./configure ", heredoc, envs)

    elif target == "upp":
        envs["WRF_ESMF_LIBS"] = envs["WRF_ESMF_LIBS_MPI"]
        envs["UPP_CONFIGURE_OPT"] = envs["UPP_CONFIGURE_MPI_OPT"]
        if any("nompi" in x for x in options):
            envs["WRF_ESMF_LIBS"] = envs["WRF_ESMF_LIBS_NOMPI"]
            envs["UPP_CONFIGURE_OPT"] = envs["UPP_CONFIGURE_NOMPI_OPT"]
        envs["JASPERLIB"] = envs["UPP_GRIB2_LIBS"]
        envs["JASPERINC"] = envs["UPP_GRIB2_INC"]
        heredoc = "<< EOF\n" + envs["UPP_CONFIGURE_OPT"] + "\nEOF"
        rc = u.run_configure("./configure ", heredoc, envs)

    elif target == "wps":
        envs["JASPERLIB"] = envs["WPS_GRIB2_LIBS"]
        envs["JASPERINC"] = envs["WPS_GRIB2_INC"]

        heredoc = "<<EOF \n" + envs["WPS_CONFIGURE_MPI_OPT"] + "\nEOF"
        rc = u.run_configure("./configure ", heredoc, envs)
        if os.environ.get("DEBUG_BUILD") is not None:
            rc = 0
        if not os.path.isfile("configure.wps"):
            rc = 1

        if any("debug" in x for x in options):
            u.sed_inplace(
                "configure.wps", r"^CFLAGS.*", "CFLAGS = " + envs["WPS_DEBUG_CFLAGS"]
            )
            u.sed_inplace(
                "configure.wps", r"^FFLAGS.*", "FFLAGS = " + envs["WPS_DEBUG_FFLAGS"]
            )
            u.sed_inplace(
                "configure.wps",
                r"^F77FLAGS.*",
                "F77FLAGS = " + envs["WPS_DEBUG_F77FLAGS"],
            )
        rc  = 0
    
    elif target == "lis":
        rc = configure_lis(my_build)
                      
    elif "wrf" in target:
        if envs["WRF_LIS"] == "1":
            rc = configure_lis(my_build)
        rc = configure_wrf(my_build)

    elif target == "utils":
        if "debug" in options:
            envs["UTILS_DEBUG"] = "1"
        rc = 0

    else:
        return 0

    if os.environ.get("DEBUG_BUILD") is not None:
        return 0
    else:
        if rc != 0:
            logger.error(target + " configure failed!")
        return rc
