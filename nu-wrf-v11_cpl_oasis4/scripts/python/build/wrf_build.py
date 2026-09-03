import logging
import os

logger = logging.getLogger(__name__)

def files_required_by_wps(my_build):
    logger.debug("Check for WRF files required by WPS.")
    adir = os.environ.get("NUWRFDIR") + "/WRF/external/"
    files = []
    if "wps" in my_build.targets and "wrf" not in my_build.targets:
        files.append(adir + "io_grib1/libio_grib1.a")
        files.append(adir + "io_grib_share/libio_grib_share.a")
        files.append(adir + "io_int/libwrfio_int.a")
        files.append(adir + "io_netcdf/libwrfio_nf.a")
        for f in files:
            if not os.path.isfile(f):
                logger.error("Cannot find file {}, will rebuild WRF".format(f))
                my_build.targets.append("wrf")
                my_build.set_targets_options(my_build.targets)
                # The following _can_ be useful but it is generaly disruptive
                # my_build.options.append('cleanfirst')
                break


def files_required_by_upp(my_build):
    logger.debug("Check for WRF files required by UPP.")
    adir = os.environ.get("NUWRFDIR") + "/WRF/"
    files = []
    if "upp" in my_build.targets and "wrf" not in my_build.targets:
        files.append(adir + "main/libwrflib.a")
        files.append(adir + "external/io_grib_share/libio_grib_share.a")
        files.append(adir + "external/io_int/libwrfio_int.a")
        files.append(adir + "external/io_netcdf/libwrfio_nf.a")
        files.append(adir + "external/esmf_time_f90/libesmf_time.a")
        files.append(adir + "external/RSL_LITE/librsl_lite.a")
        files.append(adir + "external/fftpack/fftpack5/libfftpack.a")
        files.append(adir + "frame/module_internal_header_util.o")
        files.append(adir + "frame/pack_utils.o")
        files.append(adir + "frame/module_machine.o")
        for f in files:
            if not os.path.isfile(f):
                logger.error("Cannot find file {}, will rebuild WRF".format(f))
                my_build.targets.append("wrf")
                my_build.set_targets_options(my_build.targets)
                # The following _can_ be useful but it is generaly disruptive
                # my_build.options.append('cleanfirst')
                break

