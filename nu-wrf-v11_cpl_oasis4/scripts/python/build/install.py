#!/usr/bin/env python
import os
import logging
import shared.utils as u
import build

logger = logging.getLogger("install")


def install(my_build, target):
    from shutil import copy

    nuwrf_dir = os.environ.get("NUWRFDIR")
    if target in my_build.utils_exe:
        bld_dir = "utils"
    else:
        bld_dir = my_build.target_dir[target]
    os.chdir(nuwrf_dir + "/" + bld_dir)
    if "utils" in bld_dir:
        if "geos2wrf" in target:
            installed_files = [
                "bin/createHGT.x",
                "bin/createLANDSEA.x",
                "bin/createPRESSURE.x",
                "bin/createRH.x",
                "bin/extrapIsobaric.x",
                "bin/geos2wps.x",
                "bin/merra2wrf.x",
                "bin/splitWPS.x",
                "bin/temporalInterpolation.x",
                "bin/merra2wrf.py",
            ]
        elif "sst2wrf" in target:
            installed_files = [
                "bin/sst2wrf.x",
                "bin/sst2wrf.py",
            ]
        elif "gocart2wrf" in target:
            installed_files = [
                "bin/gocart2wrf.x",
                "gocart2wrf/scripts/proc_merra2_gocart_ges_disc.sh",
            ]
        elif "casa2wrf" in target:
            installed_files = [
                "bin/Read_CO2_conc.x",
                "bin/Read_CO2_Flux.x",
                "bin/ConvertData2Netcdf.x",
                "bin/casa2wrf.x",
            ]
        elif "lisWrfDomain" in target:
            installed_files = [
                "bin/lisWrfDomain.x",
                "bin/adjust_liswrf_domain.py",
            ]
        else:
            installed_files = build.get_expected_output(my_build, target)
    else:
        if "upp" in bld_dir:
            installed_files = [
                "bin/unipost.exe",
                "bin/ndata.exe",
                "bin/copygb.exe",
                "main/libwrflib.a",
                "external/io_grib_share/libio_grib_share.a",
                "external/io_int/libwrfio_int.a",
                "external/io_netcdf/libwrfio_nf.a",
                "external/esmf_time_f90/libesmf_time.a",
                "external/RSL_LITE/librsl_lite.a",
                "external/fftpack/fftpack5/libfftpack.a",
            ]
        else:
            installed_files = build.get_expected_output(my_build, target)
    for f in installed_files:
        logger.info("Install [" + f.split("/")[-1] + "]")
        if os.path.isfile(f):
            p = os.path.dirname(my_build.prefix + "/" + bld_dir + "/" + f)
            if not os.path.exists(p):
                u.mkdir_p(p)
            logger.debug("copy " + f + " -> " + p)
            copy(f, p)
        else:
            logger.warning(f + " does not exist -> install failed.")
