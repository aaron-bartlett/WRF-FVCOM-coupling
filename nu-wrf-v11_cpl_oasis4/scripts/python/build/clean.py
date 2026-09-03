import logging
import os
from glob import glob
import shared.utils as u
import nuwrf_build as nuwrf

logger = logging.getLogger(__name__)

def clean_it(my_build, target):
    envs = my_build.env_vars

    logger.debug("Cleaning [" + target + "]")

    rc = 0

    if "allclean" in my_build.targets:
        if os.path.isfile(os.environ.get("NUWRFDIR") + "/.build_settings"):
            os.remove(os.environ.get("NUWRFDIR") + "/.build_settings")
        if os.path.isfile(os.environ.get("NUWRFDIR") + "/.modules"):
            os.remove(os.environ.get("NUWRFDIR") + "/.modules")

    if target == "ldt":
        if not os.path.isfile("touch make/configure.ldt"):
            u.sp_call("touch make/configure.ldt")
        rc = u.run_shell_command("make -C make/MAKDEP clean", envs)
        rc = u.run_shell_command("make -C make realclean", envs)
        if os.path.isfile("LDT"):
            os.remove("LDT")
        if os.path.isfile("make/configure.ldt"):
            os.remove("make/configure.ldt")
        if os.path.isfile("make/LDT_NetCDF_inc.h"):
            os.remove("make/LDT_NetCDF_inc.h")
        if os.path.isfile("make/LDT_misc.h"):
            os.remove("make/LDT_misc.h")

    if target == "lvt":
        if not os.path.isfile("touch make/configure.lvt"):
            u.sp_call("touch make/configure.lvt")
        rc = u.run_shell_command("make -C make/MAKDEP clean", envs)
        rc = u.run_shell_command("make -C make realclean", envs)
        if os.path.isfile("LVT"):
            os.remove("LVT")
        if os.path.isfile("make/configure.lvt"):
            os.remove("make/configure.lvt")
        if os.path.isfile("make/LVT_NetCDF_inc.h"):
            os.remove("make/LVT_NetCDF_inc.h")
        if os.path.isfile("make/LVT_misc.h"):
            os.remove("make/LVT_misc.h")

    if target == "doc":
        rc = u.run_shell_command("make -C tutorial clean", envs)
        rc = u.run_shell_command("make -C userguide clean", envs)

    if target in ["gsdsu", "gmp", "grad"]:
        rc = u.run_shell_command("make clean", envs)

    if target == "met":
        rc = u.run_shell_command("make clean", envs)
        rc = u.sp_call("git clean -fd . > /dev/null")

    if target == "utils" or target in nuwrf.NuwrfBuild().utils_exe:
        if target == "utils":
            rc = u.run_shell_command(
                "make allclean CONFIG_DIR=" + os.getcwd() + "/config", envs
            )
        else:
            rc = u.run_shell_command(
                "make -C shared clean CONFIG_DIR=" + os.getcwd() + "/config", envs
            )
            rc = u.run_shell_command(
                "make -C " + target + " clean CONFIG_DIR=" + os.getcwd() + "/config",
                envs,
            )

    if target in ["arw", "rip", "upp", "wps", "wrf"]:

        rc = u.run_shell_command("./clean -a", envs)

        if target == "arw":
            if os.path.isfile("configure.arwp.backup"):
                os.remove("configure.arwp.backup")

        if target == "rip":
            if os.path.isfile("configure.rip.backup"):
                os.remove("configure.rip.backup")

        if target == "upp":
            if os.path.isfile("configure.upp"):
                os.remove("configure.upp")

        if target == "wps":
            if os.path.isfile("configure.wps.backup"):
                os.remove("configure.wps.backup")

        if target == "wrf":
            if os.path.isfile("Registry/Registry.backup"):
                os.remove("Registry/Registry.backup")
            if os.path.isfile("configure.wrf.backup"):
                os.remove("configure.wrf.backup")
            for f in glob("run/namelist.input.backup.*"):
                os.remove(f)
            rc = u.run_shell_command("./clean -a", envs)

    if target == "lis":
        if os.path.isfile("LIS"):
            os.remove("LIS")
        os.chdir("make")
        u.sp_call("touch configure.lis")
        rc = u.run_shell_command("make realclean", envs)
        os.remove("configure.lis")
        os.chdir(os.environ.get("NUWRFDIR") + "/" + my_build.target_dir[target])
        if os.path.isfile(os.environ.get("NUWRFDIR")+"/.lis"):
            os.remove(os.environ.get("NUWRFDIR")+"/.lis")
        
    if "allclean" not in my_build.targets:
        my_build.remove_built_component_opt(target)

    return rc


def clean_sub(my_build, target, adir):
    t_dir = os.environ.get("NUWRFDIR") + "/" + adir
    if not os.path.isdir(t_dir):
        logger.debug("({}) does not exist.".format(t_dir))
        return
    else:
        os.chdir(t_dir)
    clean_it(my_build, target)
    os.chdir(os.environ.get("NUWRFDIR"))
