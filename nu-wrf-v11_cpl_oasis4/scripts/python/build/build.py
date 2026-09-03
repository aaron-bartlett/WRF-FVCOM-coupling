import clean
import configure
import logging
import os
import time
import shared.utils as u
import nuwrf_build as nuwrf

logger = logging.getLogger(__name__)
# For WRF ideal builds
ideal_case_names = {
    "ideal_b_wave": "em_b_wave",
    "ideal_convrad": "em_convrad",
    "ideal_heldsuarez": "em_heldsuarez",
    "ideal_les": "em_les",
    "ideal_quarter_ss": "em_quarter_ss",
    "ideal_scm_xy": "em_scm_xy",
    "ideal_scm_lis_xy": "em_scm_lis_xy",
    "ideal_tropical_cyclone": "em_tropical_cyclone",
}


def get_expected_output(my_build, target):
    envs = my_build.env_vars
    options = my_build.options
    expected_output = list()

    if target == "arw":
        expected_output.append("ARWpost.exe")

    elif target == "doc":
        expected_output.append("tutorial/tutorial.pdf")
        expected_output.append("userguide/nuwrf_userguide.pdf")

    elif target == "gmp":
        expected_output.append("../QRUN/GMP.x")

    elif target == "grad":
        expected_output.append("../QRUN/GRAD.x")
        expected_output.append("../QRUN/HOW_MANY_CPU_GRAD")

    elif target == "gsdsu":
        expected_output.append("../QRUN/GSDSU.x")
        expected_output.append("../QRUN/HOW_MANY_CPU_GSDSU")

    elif target == "ldt":
        expected_output.append("LDT")

    elif target == "lvt":
        expected_output.append("LVT")

    elif target == "met":
        expected_output.append("bin/ascii2nc")
        expected_output.append("bin/ensemble_stat")
        expected_output.append("bin/gen_vx_mask")
        expected_output.append("bin/grid_stat")
        expected_output.append("bin/gsid2mpr")
        expected_output.append("bin/gsidens2orank")
        expected_output.append('bin/lidar2nc')
        expected_output.append("bin/madis2nc")
        expected_output.append("bin/mode")
        expected_output.append("bin/mode_analysis")
        expected_output.append("bin/modis_regrid")
        expected_output.append("bin/mtd")
        expected_output.append("bin/pb2nc")
        expected_output.append("bin/pcp_combine")
        expected_output.append("bin/plot_data_plane")
        expected_output.append("bin/plot_point_obs")
        expected_output.append('bin/plot_mode_field')
        expected_output.append("bin/point_stat")
        expected_output.append("bin/regrid_data_plane")
        expected_output.append("bin/series_analysis")
        expected_output.append("bin/shift_data_plane")
        expected_output.append("bin/stat_analysis")
        expected_output.append("bin/tc_dland")
        expected_output.append("bin/tc_pairs")
        expected_output.append("bin/tc_stat")
        expected_output.append("bin/wavelet_stat")
        expected_output.append("bin/wwmca_plot")
        expected_output.append("bin/wwmca_regrid")

    elif target == "rip":
        expected_output.append("rip")
        expected_output.append("ripcomp")
        expected_output.append("ripcut")
        expected_output.append("ripdp_mm5")
        expected_output.append("ripdp_wrfarw")
        expected_output.append("ripdp_wrfnmm")
        expected_output.append("ripinterp")
        expected_output.append("ripshow")
        expected_output.append("showtraj")
        expected_output.append("tabdiag")
        expected_output.append("upscale")

    elif target == "upp":
        expected_output.append("bin/unipost.exe")
        expected_output.append("bin/ndate.exe")
        expected_output.append("bin/copygb.exe")

    elif target == "utils":
        for i in nuwrf.NuwrfBuild().utils_exe:
            expected_output.append("bin/" + i + ".x")

    elif target in nuwrf.NuwrfBuild().utils_exe:
        expected_output.append("bin/" + target + ".x")

    elif target == "wps":
        # Core WPS
        expected_output.append("geogrid.exe")
        expected_output.append("metgrid.exe")
        expected_output.append("ungrib.exe")
        # WPS Util
        expected_output.append("util/avg_tsfc.exe")
        expected_output.append("util/calc_ecmwf_p.exe")
        expected_output.append("ungrib/g1print.exe")
        expected_output.append("ungrib/g2print.exe")
        expected_output.append("util/height_ukmo.exe")
        expected_output.append("util/int2nc.exe")
        expected_output.append("util/mod_levs.exe")
        expected_output.append("util/plotfmt.exe")
        expected_output.append("util/plotgrids.exe")
        expected_output.append("util/rd_intermediate.exe")

    elif target == "lis":
        expected_output.append("LIS")

    elif target == "wrf":
        # Specify directories for each ideal case target
        if "ideal_case" in options:
            expected_output.append("main/ideal.exe")
            expected_output.append("main/wrf.exe")
        else:
            expected_output.append("main/ndown.exe")
            expected_output.append("main/real.exe")
            expected_output.append("main/tc.exe")
            expected_output.append("main/wrf.exe")
            if envs["WRF_CHEM"] == "1":
                expected_output.append("chem/convert_emiss.exe")

    return expected_output


def is_built(target, expected_output):
    all_built = len(expected_output)
    i = 0
    if i < all_built:
        for executable in expected_output:
            if not os.path.isfile(executable):
                break
            else:
                i = i + 1
                continue
        if all_built == i:
            logger.info(target + " is already built.")
            return True
    return False


def build_lis(my_build, t):
    logger.info(" - building [" + t + "] -")
    make_dir = os.path.join(os.environ.get("NUWRFDIR"),
                            my_build.target_dir["lis"], "make")
    os.chdir(make_dir)
    rc = u.sp_check_call_make_log("make -f Makefile "+t+" -j 4", my_build.env_vars)
    if rc == 0 :
        u.sp_call("touch "+os.environ.get("NUWRFDIR")+"/.lis")
    
def build_wrf(my_build):
    os.chdir(os.environ.get("NUWRFDIR") + "/" + my_build.target_dir["wrf"])

    Jm = my_build.env_vars["J"]
    if my_build.ideal_case:
        idn = my_build.env_vars["ideal_case_name"]
        # Remove all symbolic links to ideal.exe
        fs = u.find_files("test", "\*.exe")
        for f in fs:
            if os.environ.get("DEBUG_BUILD") is None:
                os.remove(f)
        logger.info(" - building ideal case [" + idn + "] -")
        u.sp_check_call_make_log(
            "./compile " + Jm + " " + ideal_case_names[idn], my_build.env_vars
        )
    else:
        logger.info(" - building [em_real] - ")
        if my_build.env_vars["WRF_CHEM"] == "0":
            u.sp_check_call_make_log("./compile "+Jm+" em_real", my_build.env_vars)
        else:
            logger.info(" --- building chem support -")
            u.sp_check_call_make_log("./compile "+Jm+" em_real", my_build.env_vars)
            u.sp_check_call_make_log("./compile "+Jm+" emi_conv", my_build.env_vars)


def compile_it(my_build, target):
    logger.info("Compile [" + target + "]")

    expected_output = get_expected_output(my_build, target)    
    if "rebuild" in my_build.options:  # and 'wrf' in target:
        for executable in expected_output:
            if os.environ.get("DEBUG_BUILD") is None:
                try:
                    os.remove(executable)
                except OSError:
                    pass
    # Else, check if we are all done and return to calling routine
    else:


        
        if "cleanfirst" in my_build.options:
            for f in expected_output:
                if os.environ.get("DEBUG_BUILD") is None:
                    try:
                        os.remove(f)
                    except OSError:
                        pass
        else:
            if is_built(target, expected_output):
                return
    envs = my_build.env_vars

    if target in ["arw", "doc", "rip", "upp"]:
        u.sp_check_call_make_log("./compile", envs)

    elif target in ["ldt", "lvt"]:
        u.sp_check_call_make_log("./compile -j 4", envs)
        
    elif target in ["gmp", "grad", "gsdsu"]:
        STR = target.upper()
        make_args = (
            "make -f "
            + envs[STR+"_MAKEFILE"]
            + " INC_NETCDF="
            + envs[STR+"_NETCDF_INCDIR"]
            + " LD_NETCDF="
            + envs[STR+"_NETCDF_LIBDIR"]
            + " NETCDF4_DEP_LIB="
            + "'"
            + envs[STR+"_NETCDF4_DEP_LIB"]
            + "'"
            + " all"
        )
        u.sp_check_call_make_log(make_args, envs)

    elif target == "met":
        u.sp_check_call_make_log("make install", envs)

    elif target == "wps":
        u.sp_check_call_make_log("./compile wps", envs)
        u.sp_check_call_make_log("./compile util", envs)

    elif target == "lis":
        my_build.env_vars["BUILD_LIS"] = "1"
        build_lis(my_build, "LIS")
    
    elif target == "wrf":
        # This is used is LIS's Makefile to avoid rebuilding LIS code
        my_build.env_vars["BUILD_LIS"] = "0"
        if "wrfonly" not in my_build.targets and not os.path.exists(os.environ.get("NUWRFDIR")+"/.lis"):
            my_build.env_vars["BUILD_LIS"] = "1"
        # CPP flag in WRF code to build with LIS coupling
        if my_build.env_vars["WRF_LIS"] == "1":
            build_lis(my_build, "explis")
        build_wrf(my_build)

    if os.environ.get("DEBUG_BUILD") is not None:
        logger.info("[" + target + "] build was successful.")
        return

    # Assume nothing has built
    all_built = False
    num_built = len(expected_output)
    count = 0
    for f in expected_output:
        if not os.path.isfile(f):
            logger.warning("--- Target: {} was not built!".format(f))
        else:
            count = count + 1
        if count == num_built:
            all_built = True
    if not all_built:
        open(os.environ.get("NUWRFDIR") + "/.build_failed", "a").close()
        logger.error("[" + target + "] build FAILED")
        return

    logger.info("[" + target + "] build was successful.")


def compile_utils(my_build, target):
    logger.info("Compile [" + target + "]")
    envs = my_build.env_vars
    options = my_build.options

    make_target = target
    if "utils" in target:
        make_target = "all"

    J = envs["J"]
    if "prep_chem_sources" in make_target or "all" in make_target:
        J = ""

    # Required files:
    expected_output = []
    if "all" in make_target:
        for f in nuwrf.NuwrfBuild().utils_exe:
            expected_output.append("bin/" + f + ".x")
    else:
        expected_output.append("bin/" + make_target + ".x")

    if "rebuild" in my_build.options:  # and 'wrf' in target:
        for executable in expected_output:
            if os.environ.get("DEBUG_BUILD") is None:
                try:
                    os.remove(executable)
                except OSError:
                    pass
    # Else, check if we are all done and return to calling routine
    else:
        if "cleanfirst" in my_build.options:
            for f in expected_output:
                if os.environ.get("DEBUG_BUILD") is None:
                    try:
                        os.remove(f)
                    except OSError:
                        pass
        else:
            if is_built(target, expected_output):
                return

    if "debug" in options:
        envs["UTILS_DEBUG"] = "yes"
    make_args = (
        "make " + J + " " + make_target + " CONFIG_DIR=" + os.getcwd() + "/config"
    )
    rc = u.sp_check_call_make_log(make_args, envs)

    if os.environ.get("DEBUG_BUILD") is not None:
        logger.info("[" + target + "] build was successful.")
        return
    for f in expected_output:
        if not os.path.exists(os.path.join(os.getcwd(),f)):
            logger.error("[" + target + "] build FAILED: {} was not built!".format(f))
            return
    else:
        if rc != 0:
            open(os.environ.get("NUWRFDIR") + "/.build_failed", "a").close()
            logger.error("utils build failed!")
            return

    logger.info("[" + target + "] build was successful.")


def main(my_build, target):
    start_time = time.time()
    nuwrf_dir = os.environ.get("NUWRFDIR")

    if target in my_build.utils_exe:
        build_dir = nuwrf_dir + "/utils"
    else:
        build_dir = nuwrf_dir + "/" + my_build.target_dir[target]

    os.chdir(build_dir)
    logger.debug("Entering " + os.getcwd())

    build_options = my_build.options
    if "cleanfirst" in build_options:
        rc = clean.clean_it(my_build, target)
        if os.environ.get("DEBUG_BUILD") is None:
            if rc != 0:
                logger.error("Clean failed")
                return

    rc = configure.configure(my_build, target)
    if os.environ.get("DEBUG_BUILD") is None:
        if rc != 0:
            logger.error("Configure failed")
            return
    
    if target == "utils" or target in my_build.utils_exe:
        os.chdir(nuwrf_dir + "/utils")
        compile_utils(my_build, target)
    else:
        compile_it(my_build, target)

    # Done
    end_time = time.time() - start_time
    logger.info("[" + target + "] build time taken = %f" % end_time)
