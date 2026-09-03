"""
Specific functions for NU-WRF regression testing.
"""
import os
import stat
import glob
import re
import time
import subprocess as sp
import logging
from shutil import copytree, ignore_patterns

try:
    import ConfigParser
except ImportError:
    from configparser import ConfigParser
import platform
import sys
import shared.utils as utils
import RegTest
import RegPool as pool

logger = logging.getLogger("tools")

"""
NU-WRF regression testing is divided into three stages:

1) Initialization of the testing environment
   reg_initialize
       create_test_list
       setup_testing_environment
           git clone or svn checkout repository
           create_test_dir_structure
2) Running all the tests
   reg_run
       reg_builds
           create_build_commands
               create_build_script(s)
               run_commands concurrently
       reg_runs
           create_run_commands
               create_run_script(s)
               run_commands concurrently
3) Creating a test report

"""


def reg_initialize(config, comp_config, time_stamp):
    task_list = create_test_list(config, time_stamp)
    setup_testing_environment(config, comp_config, task_list)
    return task_list


def reg_run(config, comp_config, task_list, time_stamp):
    build_results = reg_builds(config, comp_config, task_list, time_stamp)
    reg_runs(config, comp_config, task_list, build_results)
    return None


def reg_finalize(config, comp_config, end_time):
    user_config = utils.config_section_map(config, "USERCONFIG")
    # Notify results via email
    if user_config["mail_to"]:
        send_test_report(config, comp_config, end_time)
    return None


# --- INITIALIZE ---
def setup_testing_environment(config, comp_config, task_list):
    logger.info("Setup testing environment...")
    user_config = utils.config_section_map(config, "USERCONFIG")

    if "yes" in user_config["clean_scratch"]:
        utils.clean_scratch(config)

    if not user_config["model_dir"]:
        if user_config["repo_type"] == "git":
            git_clone_repository(config)
        else:
            svn_checkout_repository(config)
    create_test_dir_structure(config, comp_config, task_list)


def create_test_list(config, time_stamp):
    """
    Create/return a list of model run configurations specified in config file
    The run configurations are a list of RegRun instances.
    """
    sections = config.sections()
    # Store model configurations in a list and let each item in the list
    # have access to the user-defined options.
    test_list = list()
    # NOTE: Build tasks are the first element in the task list
    build_tasks = set()
    # Retrieve all model run configurations. For convenience divide the sections
    # in the configuration file into two types: CONFIG and others. The former
    # have CONFIG in their names. Thus, if a section name does NOT have CONFIG
    # in its name then it is a model configuration (or test case).
    for test_case in sections:
        match = not re.search("CONFIG", test_case)
        build = test_case.split("_", 1)[0]
        if "wrflis" in test_case:
            if "scm" in test_case:
                build = "scm"
            else:
                build = "wrf"
        if match:
            for compiler in config.get(test_case, "compilers").split(","):
                # Create a build task for each "compiler" in "compilers"
                build_tasks.add((compiler, build))
                # Get test_case settings, in a ConfigParser object
                reg = get_reg_config(config, test_case)
                # Create a regression test instance
                rt = RegTest.RegTest(test_case)
                # Set its compiler
                reg.set(test_case, "compiler", compiler)
                # Each test case has the same time stamp
                reg.set(test_case, "time_stamp", time_stamp)
                # Use ConfigParser object, reg, to populate test instance data
                rt.set_opts(reg, test_case)
                # Finally, append to test_list
                test_list.append(rt)

    # Note that build tasks are the first element of the task list
    test_list.insert(0, build_tasks)
    return test_list


def get_reg_config(config, test_name):
    # From the config file, get name, values from the 'test_name' section
    # and return a ConfigParser with those settings.
    try:
        reg_config = ConfigParser.ConfigParser()
    except AttributeError:
        reg_config = ConfigParser()
    reg_config.add_section(test_name)
    for name, value in config.items(test_name):
        reg_config.set(test_name, name, value)
    for name, value in config.items("USERCONFIG"):
        reg_config.set(test_name, name, value)
    return reg_config


def svn_checkout_repository(config):
    logger.info("Checkout user-specified svn repository...")
    user_config = utils.config_section_map(config, "USERCONFIG")
    repo = user_config["repo_url"]

    svn_repo = user_config["scratch_dir"] + "/builds/nu-wrf_repo"
    if os.path.exists(svn_repo):
        return
    cmd = ["svn", "--quiet", "checkout", repo, svn_repo]

    logger.info("Checking out %s into %s", repo, svn_repo)
    cwd = os.getcwd()
    proc = sp.Popen(cmd)
    proc.wait()
    os.chdir(cwd)


def git_clone_repository(config):
    logger.info("Clone user-specified git repository...")
    user_config = utils.config_section_map(config, "USERCONFIG")
    temp = user_config["scratch_dir"]
    repo = user_config["repo_url"]
    branch = user_config["repo_branch"]

    clone = temp + "/builds/nu-wrf_repo"
    if not branch:
        cmd = ["git", "clone", "-q", repo, clone]
    else:
        cmd = ["git", "clone", "-q", "-b", branch, repo, clone]

    results_dir = temp + "/results/"
    if os.path.exists(clone):
        os.chdir(clone)
        cmd = "git log --pretty=format:'%h - %an, %ar : %s' --since=1.day"
        os.system(cmd + ">" + results_dir + "gitLog")
        return

    logger.info("Cloning %s into %s", repo, clone)
    cwd = os.getcwd()
    proc = sp.Popen(cmd)
    proc.wait()
    os.chdir(clone)
    cmd = "git log --pretty=format:'%h - %an, %ar : %s' --since=1.day"
    os.system(cmd + ">" + results_dir + "gitLog")
    os.chdir(cwd)


def create_test_dir_structure(config, comp_config, test_list):
    logger.info("Create directories for each NU-WRF build...")
    user_config = utils.config_section_map(config, "USERCONFIG")
    builds_dir = user_config["scratch_dir"] + "/builds"
    results_dir = user_config["scratch_dir"] + "/results"

    time_stamp = test_list[1].get_opt("time_stamp")

    build_tasks = test_list[0]
    for item in build_tasks:
        compiler = item[0]
        build = item[1]
        logger.debug("Create " + build + " build task using " + compiler + " settings")
        available = get_compilers(comp_config)
        # Sanity check
        if compiler.split("-")[0] in available:
            build_dir = builds_dir + "/" + compiler + "." + time_stamp + "/" + build
            if user_config["model_dir"]:
                build_dir_link = builds_dir + "/" + compiler + "." + time_stamp + "/"
                utils.mkdir_p(build_dir_link)
                repo_dir = user_config["model_dir"]
                logger.debug("Creating symlink " + repo_dir + " -> " + build_dir)
                os.symlink(repo_dir, build_dir)
            else:
                logger.debug("Creating build dir " + build_dir)
                repo_dir = builds_dir + "/nu-wrf_repo"
                copytree(
                    repo_dir,
                    build_dir,
                    symlinks=True,
                    ignore=ignore_patterns(
                        "GSDSU",
                        "MET",
                        "UPP",
                        "docs",
                        "." + user_config["repo_type"],
                    ),
                )
        else:
            logger.error("Unable to create directory for " + build + " build")

    # Create  directories to hold NU-WRF results for each run
    iter_runs = iter(test_list)
    next(iter_runs)
    for t in iter_runs:
        compiler = t.get_opt("compiler")
        run_dir = results_dir + "/" + compiler + "." + time_stamp + "/" + t.name
        logger.debug("Creating results dir " + run_dir)
        utils.mkdir_p(run_dir)


# --- RUN ---
def reg_builds(config, comp_config, test_list, time_stamp):
    user_config = utils.config_section_map(config, "USERCONFIG")
    builds_dir = user_config["scratch_dir"] + "/builds"

    # Accumulate build results
    build_results = dict()
    if not user_config["model_dir"]:
        build_commands = create_build_commands(config, comp_config, test_list)
        if build_commands:
            logger.info("Compiling...")
            if "DEBUG_REG" in os.environ:
                logger.info(build_commands)
            else:
                if "Darwin" in platform.system():
                    pool.run_commands(build_commands, "no")
                else:
                    pool.run_commands(build_commands, user_config["use_batch"])

        for item in test_list[0]:
            compiler = item[0]
            build = item[1]
            available = get_compilers(comp_config)
            report_file = (
                builds_dir + "/" + compiler + "." + time_stamp + "/builds.diff"
            )
            if "DEBUG_REG" in os.environ:
                build_results[(build, compiler)] = "b+"
                logger.info("Build " + build + " was successful.")
                results = [
                    build,
                    compiler.center(10, " "),
                    "b+".center(10, " "),
                    "-".center(10, " "),
                ]
            else:
                # Sanity check
                if compiler.split("-")[0] in available:
                    build_dir = (
                        builds_dir + "/" + compiler + "." + time_stamp + "/" + build
                    )
                    results = [
                        build,
                        compiler.center(10, " "),
                        "-".center(10, " "),
                        "-".center(10, " "),
                    ]
                    if os.path.exists(build_dir + "/.build_ok"):
                        build_results[(build, compiler)] = "b+"
                        logger.info("Build " + build + " was successful.")
                        results = [
                            build,
                            compiler.center(10, " "),
                            "b+".center(10, " "),
                            "-".center(10, " "),
                        ]
                    else:
                        build_results[(build, compiler)] = "bF"
                        logger.info("Build " + build + " FAILED.")
                        results = [
                            build,
                            compiler.center(10, " "),
                            "bF".center(10, " "),
                            "-".center(10, " "),
                        ]
            with open(report_file, "a") as fp:
                utils.write_results(results, fp)

    return build_results


def create_build_commands(config, comp_config, test_list):
    user_config = utils.config_section_map(config, "USERCONFIG")
    builds_dir = user_config["scratch_dir"] + "/builds"

    time_stamp = test_list[1].get_opt("time_stamp")

    build_commands = []
    # Create build commands for each NU-WRF build
    for item in test_list[0]:
        comp_mpi = item[0]
        build = item[1]
        available = get_compilers(comp_config)
        if comp_mpi.split("-")[0] in available:
            build_dir = builds_dir + "/" + comp_mpi + "." + time_stamp + "/" + build
            ideal_build = True if "scm" in build else False
            script_command = create_build_script(
                config, comp_config, build, build_dir, ideal_build, comp_mpi
            )
            if script_command:
                build_commands.append(script_command)
    return build_commands


def create_build_script(config, comp_config, build, build_dir, ideal_build, comp_mpi):
    """ Create a build script file and return a command to run it."""
    user_config = utils.config_section_map(config, "USERCONFIG")
    filename = build_dir + "/" + build + ".bash"
    logger.debug("Creating build script in: " + build_dir)
    with open(filename, "w") as fp:
        fp.write("#!/bin/bash" + "\n")
        add_batch_directives(config, fp, build, build_dir, is_type=0)
        fp.write("cd " + build_dir + "\n")
        build_arg = "lis,wrf,wps,ldt,arw,rip,geos2wrf,sst2wrf,lisWrfDomain"
        if ideal_build:  # Only one ideal case is tested
            build_arg = "-o nompi lis,ideal_scm_lis_xy,wps,ldt,lis4scm"
        else:
            if "chem" in build:
                build_arg = "lis,chem,wps,ldt,utils"
            elif "kpp" in build:
                build_arg = "lis,kpp,wps,prep_chem_sources,gocart2wrf"
        logger.debug(build + " build components: " + build_arg)

        comp_config_section = utils.config_section_map(comp_config, "COMPCONFIG")
        comp_mpi_vendor = "intel-intelmpi"
        for mod in comp_config_section["modulelist"].split(","):
            if comp_mpi in mod:
                comp_mpi_vendor = mod

        libdir_tag = os.path.join(user_config["nuwrflib_dir"], comp_mpi_vendor)
        fp.write("export LIBDIR_TAG=" + libdir_tag + "\n")
        fp.write("export COMPILER_VENDOR=" + comp_mpi_vendor.split("-")[0] + "\n")
        fp.write("export MPI_VENDOR=" + comp_mpi_vendor.split("-")[1] + "\n")

        fp.write("./build.sh " + build_arg + "\n")
        fp.write("if [ ! -f ./.build_failed ]; then touch .build_ok; exit; fi" + "\n")

	# *** Temporary DISCOVER SLES12 hack ***
        # Try to recompile again -up to 4 times- then give up.
        if "intel" in comp_mpi_vendor and "discover" in platform.node():
            fp.write("rm -f ./.build_failed" + "\n")
            fp.write("counter=0" + "\n")
            fp.write("while  [ $counter -lt 3 ]; do" + "\n")
            fp.write("./build.sh " + build_arg + "\n")
            fp.write("if [ ! -f ./.build_failed ]; then touch .build_ok; exit; fi" + "\n")
            fp.write("let counter=counter+1" + "\n")
            fp.write("done" + "\n")

    if "yes" in user_config["use_batch"]:
        script_command = "cd " + build_dir + ";sbatch " + build + ".bash"
    else:
        if "Darwin" in platform.system():
            script_command = "chmod +x " + filename + ";" + filename
        else:
            script_command = "chmod +x " + filename + ";ssh discover-sles12 " + filename

    return script_command


def reg_runs(config, comp_config, test_list, build_results):
    user_config = utils.config_section_map(config, "USERCONFIG")
    run_commands = create_run_commands(config, comp_config, test_list, build_results)
    try:
        if "yes" in user_config["setup_runs"]:
            logger.info("Setup is DONE.")
            sys.exit()
    except KeyError as e:
        logger.error(str(e) + " option does not exist in your config file.")
        sys.exit()
    if run_commands:
        logger.info("Running test cases...")
        if "DEBUG_REG" in os.environ:
            logger.info(run_commands)
        else:
            pool.run_commands(run_commands, user_config["use_batch"])


def create_run_commands(config, comp_config, test_list, build_results):
    user_config = utils.config_section_map(config, "USERCONFIG")
    run_commands = []

    # Skip first item which contains build tasks.
    iter_runs = iter(test_list)
    next(iter_runs)
    for t in iter_runs:
        if "compile_only" in t.get_opt("verification"):
            continue
        build = t.name.split("_", 1)[0]
        if "wrflis" in t.name:
            if "scm" in t.name:
                build = "scm"
            else:
                build = "wrf"
        run_command = None
        if user_config["model_dir"]:
            compilers = t.get_opt("compiler").split(",")
            for compiler in compilers:
                run_command = create_run_script(config, comp_config, t, compiler)
        else:
            for build_set, result in build_results.items():
                build_ok = "b+" in result
                is_the_build = build == build_set[0]
                is_the_compiler = t.get_opt("compiler") in build_set
                built = build_ok and is_the_build and is_the_compiler
                if built or "DEBUG_REG" in os.environ:
                    run_command = create_run_script(
                        config, comp_config, t, build_set[1]
                    )

        if run_command:
            run_commands.append(run_command)

    return run_commands


def create_run_script(config, comp_config, run, comp):
    """ This function creates a run script file and returns a command to executed."""
    user_config = utils.config_section_map(config, "USERCONFIG")
    logger.debug("Creating run script for test: " + run.name)
    time_stamp = run.get_opt("time_stamp")
    results_dir = run.get_opt("scratch_dir") + "/results/" + comp + "." + time_stamp
    run_dir = results_dir + "/" + run.name
    build_dir = run.get_opt("scratch_dir") + "/builds/" + comp + "." + time_stamp

    test_case = run.name.split("_", 1)
    build = test_case[0]
    if "wrflis_" in run.name:
        if "scm" in run.name:
            build = "scm"
        else:
            build = "wrf"
    if "wrf_" in run.name:
        build = "wrf"
    if "chem_" in run.name or "kpp_" in run.name:
        if "chem_" in run.name:
            build = "chem"
        else:
            build = "kpp"
    model_dir = run.get_opt("model_dir")
    if model_dir:
        build_dir = model_dir
    else:
        build_dir = build_dir + "/" + build
    exe_dir = None
    if user_config["exe_dir"]:
        exe_dir = user_config["exe_dir"]

    data_dir = run.get_opt("data_dir") + "/" + "/".join(test_case)
    # This is hardwired relative to NUWRFDIR
    test_dir = build_dir + "/testcases/" + "/".join(test_case)
    filename = results_dir + "/" + run.name + "/" + run.name + ".bash"
    with open(filename, "w") as fp:

        fp.write("#!/bin/bash" + "\n")
        npes_list = [int(x) for x in run.get_opt("npes")]
        add_batch_directives(
            config,
            fp,
            test_case[1],
            run_dir,
            npes=max(npes_list),
            is_type=1,
            wall_time=run.get_opt("wall_time"),
            qos=run.get_opt("qos"),
        )
        fp.write("cd " + run_dir + "\n")
        fp.write("# source common file for modules and paths \n")
        fp.write("source ./common.reg || exit 1 \n")
        fp.write("if [ -d " + data_dir + " ]; then \n")
        fp.write("ln -sf " + data_dir + "\n")
        fp.write("fi \n")
        fp.write(
            "python "
            + run.get_opt("scripts_dir")
            + "/"
            + "regression.py "
            + run.name
            + "\n"
        )

        fp.write(" " + "\n")
        fp.close()

    create_common_script(
        config, comp_config, comp, run_dir, build_dir, data_dir, exe_dir
    )
    create_component_script(
        config, comp_config, run, comp, test_dir, run_dir, build_dir, data_dir, exe_dir
    )
    create_run_config(results_dir, run)

    if "yes" in run.get_opt("use_batch"):
        run_command = "sbatch " + filename
    else:
        run_command = "chmod +x " + filename + "; " + filename

    return run_command


def create_common_script(
    config, comp_config, comp, run_dir, build_dir, data_dir, exe_dir=None
):
    user_config = utils.config_section_map(config, "USERCONFIG")
    comp_config_section = utils.config_section_map(comp_config, "COMPCONFIG")
    with open(run_dir + "/common.reg", "w") as fp:
        # Module environment
        comp_mpi_vendor = "intel-intelmpi"
        if "yes" in user_config["use_modules"]:
            if "Darwin" in platform.system():
                # CC: This is not portable...just my MAC so far
                fp.write(". /opt/local/share/modules/bin" + "\n")
                fp.write("module purge" + "\n")
            else:
                fp.write("umask 022" + "\n")
                fp.write(". /usr/share/modules/init/bash" + "\n")
                fp.write("module purge" + "\n")
                fp.write("unset LD_LIBRARY_PATH" + "\n")
            for mod in comp_config_section["modulelist"].split(","):
                if re.search(comp, mod):
                    comp_mpi_vendor = mod
                    for mm in comp_config_section[mod].split(","):
                        cmd = "module load " + mm + "\n"
                        fp.write(cmd)


        libdir_tag = os.path.join(user_config["nuwrflib_dir"], comp_mpi_vendor)
        fp.write("LIBDIR_TAG=" + libdir_tag + "\n")
        # Need a better way to do this
        if "Linux" in platform.system():
            fp.write("export LD_LIBRARY_PATH=$LIBDIR_TAG/jasper/lib:$LD_LIBRARY_PATH" + "\n")
            fp.write("export LD_LIBRARY_PATH=$LIBDIR_TAG/gdal/lib:$LD_LIBRARY_PATH" + "\n")
            fp.write("export LD_LIBRARY_PATH=$LIBDIR_TAG/jpeg/lib:$LD_LIBRARY_PATH" + "\n")
            fp.write("export LD_LIBRARY_PATH=$LIBDIR_TAG/hdf5/lib:$LD_LIBRARY_PATH" + "\n")
            fp.write("export LD_LIBRARY_PATH=$LIBDIR_TAG/netcdf4/lib:$LD_LIBRARY_PATH" + "\n")
            fp.write("export LD_LIBRARY_PATH=$LIBDIR_TAG/eccodes/lib:$LD_LIBRARY_PATH" + "\n")
        if "Linux" in platform.system():
            fp.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64" + "\n")
        fp.write("LISDIR=" + user_config["lis_project_dir"] + "\n")
        fp.write("DATADIR=" + data_dir + "\n")
        fp.write("RUNDIR=" + run_dir + "\n")
        fp.write("NUWRFDIR=" + build_dir + "\n")
        if exe_dir:
            fp.write("EXEDIR=" + exe_dir + "\n")
        else:
            fp.write("EXEDIR=" + build_dir + "\n")
        fp.write("export PATH=$LIBDIR_TAG/netcdf4/bin:$PATH"+"\n")


def create_component_script(
    config, comp_config, run, comp, test_dir, run_dir, build_dir, data_dir, exe_dir=None
):
    user_config = utils.config_section_map(config, "USERCONFIG")

    all_files = glob.glob(test_dir + "/*")
    [utils.copy_file(f, run_dir) for f in all_files]

    comp_index = 0
    for f in run.components.split(","):
        # Get specified NPES values used to run this component
        npes = run.get_npes_val(comp_index)

        reg_file = run_dir + "/" + f + ".reg"
        with open(reg_file, "w") as fp:
            fp.write("#!/bin/bash" + "\n")
            add_batch_directives(config, fp, f, run_dir, npes=npes, is_type=2)
            fp.write("# source common file for modules and paths \n")
            fp.write("source ./common.reg || exit 1 \n")

            if "sgimpt" in comp:
                mpicmd = "mpiexec_mpt -n " + npes
            elif "gnu" in comp or "intelmpi" in comp:
                mpicmd = "mpirun -np " + npes
            else:
                mpicmd = "mpirun -np " + npes
            fp.write("export MPIRUN='{}' \n".format(mpicmd))

            # Write script segment stored in testcases directory
            with open(test_dir + "/" + f + ".reg", "r") as ft:
                fp.write(ft.read())
        comp_index += 1
        st = os.stat(reg_file)
        os.chmod(reg_file, st.st_mode | stat.S_IEXEC)


def add_batch_directives(
    config, fp, build, run_dir, npes=8, is_type=0, wall_time="2:00:00", qos="allnccs"
):
    user_config = utils.config_section_map(config, "USERCONFIG")

    if "yes" in user_config["use_batch"]:

        sponsor_id = user_config["sponsor_id"]
        out_name = run_dir + "/" + build + ".out"
        err_name = run_dir + "/" + build + ".err"
        str_npes = str(npes)
        int_npes = int(npes)

        # Batch directives are only used if
        # 0) use_batch='yes' AND
        #    1) it's a build script or (is_type == 0)
        #    2) it's a run script or (is_type == 1)
        #    3) it's a component script (is_type == 2) AND setup_runs ='yes'
        if "discover" in platform.node():
            case_3 = is_type == 2 and "yes" in user_config["setup_runs"]
            if is_type == 0 or is_type == 1 or case_3:
                fp.write("#SBATCH -J " + build + "\n")
                fp.write("#SBATCH -o " + out_name + "\n")
                fp.write("#SBATCH -e " + err_name + "\n")
                fp.write("#SBATCH --account=" + sponsor_id + "\n")
                fp.write("#SBATCH --ntasks=" + str_npes + "\n")
                fp.write("#SBATCH --qos=" + qos + "\n")
                if 'debug' in qos:
                    fp.write("#SBATCH --time=1:00:00" + "\n")
                else:
                    fp.write("#SBATCH --time=" + wall_time + "\n")

        else:  # Must be Pleiades
            if int_npes > 24:
                if int_npes % 24 == 0:
                    nodes = int_npes / 24
                else:
                    nodes = (int_npes / 24) + 1
            else:
                nodes = 1
            fp.write(
                "#PBS -lselect="
                + str(nodes)
                + ":ncpus="
                + str_npes
                + ":mpiprocs="
                + str_npes
                + ":model=has\n"
            )
            fp.write("#PBS -q devel\n")
            fp.write("#PBS -o " + out_name + "\n")
            fp.write("#PBS -e " + err_name + "\n")
            fp.write("#PBS -W group_list=" + sponsor_id + "\n")
            wall_time = "02:00:00"
            if "scm" in build:
                wall_time = "00:10:00"
            elif "arw_" in build:
                wall_time = "00:30:00"
            fp.write("#PBS -l walltime=" + wall_time + "\n")


def create_run_config(results_dir, run):
    # Create a config file for regression.py script.
    # Note: there is one config file for each run/compiler combination
    run_dir = results_dir + "/" + run.name

    reg_config = run.test_to_parser()
    time_stamp = run.get_opt("time_stamp")
    reg_config.set(run.name, "time_stamp", time_stamp)

    reg_config.set(run.name, "components", run.get_opt("components"))

    filename = run_dir + "/" + run.name + ".cfg"
    with open(filename, "w") as fp:
        reg_config.write(fp)


# --- FINALIZE ---
def send_test_report(config, compconfig, end_time):
    # Create a report and notify via email
    user_config = utils.config_section_map(config, "USERCONFIG")
    mail_to = user_config["mail_to"]
    branch = user_config["repo_branch"]
    if not branch:
        branch = "none specified"
    results_dir = user_config["scratch_dir"] + "/results/"
    builds_dir = user_config["scratch_dir"] + "/builds/"
    build_type = user_config["build_type"]
    message = user_config["message"]
    compilers = get_compilers(compconfig)
    use_html = user_config["use_html"]
    test_report = results_dir + "/" + "test_report.txt"

    with open(test_report, "w") as fp:
        if "yes" in use_html:
            fp.write("<html><pre>\n")
        if "DEBUG_REG" in os.environ:
            fp.write(message + " - DRY RUN ONLY! \n")
        else:
            fp.write(message + " \n")
        if user_config["model_dir"]:
            fp.write("NU-WRF build directory: " + user_config["model_dir"] + "\n")
            fp.write("-" * 80 + "\n")
        else:
            fp.write("Repository: " + user_config["repo_url"] + "\n")
            fp.write("-" * 80 + "\n")
            fp.write("Branch: " + branch)
            fp.write("  --  Build type: " + build_type + "\n")
            fp.write("-" * 80 + "\n")
        fp.write(
            "{0:^40}{1:^20}{2:^10}{3:^10}\n".format(
                "BUILD and/or RUN NAME".center(40, " "),
                "COMPILER".center(20, " "),
                "RESULT".center(10, " "),
                "BASELINE",
            ).center(10, " ")
        )
        fp.write("-" * 80 + "\n")

        sp.call(
            "find "
            + builds_dir
            + " -name \*.diff -exec cat {} \; >"
            + builds_dir
            + "/"
            + "alldiffs",
            shell=True,
        )
        sp.call(
            "cat "
            + builds_dir
            + "/"
            + "alldiffs | sort -k 1,1 >"
            + builds_dir
            + "/"
            + "sorteddiffs",
            shell=True,
        )
        with open(builds_dir + "/" + "sorteddiffs", "r") as inf:
            fp.write(inf.read())
        os.remove(builds_dir + "/" + "alldiffs")
        os.remove(builds_dir + "/" + "sorteddiffs")

        sp.call(
            "find "
            + results_dir
            + " -name \*.diff -exec cat {} \; >"
            + results_dir
            + "/"
            + "alldiffs",
            shell=True,
        )
        sp.call(
            "cat "
            + results_dir
            + "/"
            + "alldiffs | sort -k 1,1 >"
            + results_dir
            + "/"
            + "sorteddiffs",
            shell=True,
        )
        with open(results_dir + "/" + "sorteddiffs", "r") as inf:
            fp.write(inf.read())
        os.remove(results_dir + "/" + "alldiffs")
        os.remove(results_dir + "/" + "sorteddiffs")

        fp.write("-" * 80 + "\n")
        hhmmss = time.strftime("%H:%M:%S", time.gmtime(end_time))
        fp.write("Time taken = %s \n" % hhmmss)
        fp.write("-" * 80 + "\n")
        fp.write("Legend:\n")
        fp.write("-" * 7 + "\n")
        fp.write("+   : task success\n")
        fp.write("F   : task failure\n")
        fp.write("C   : created baseline\n")
        fp.write("b   : _build_ task\n")
        fp.write("r   : _run_ task\n")
        fp.write("v   : _verification_ task\n")
        fp.write("-   : Not available\n")
        fp.write("Notes:\n")
        fp.write("-" * 6 + "\n")
        comp_config = utils.config_section_map(compconfig, "COMPCONFIG")
        comp_versions = comp_config["compiler_versions"].split(",")
        i = 0
        for comp in compilers:
            fp.write(comp + " compiler version: " + comp_versions[i] + "\n")
            i += 1
        fp.write("Results in: " + results_dir + "\n")
        fp.write("-" * 80 + "\n")
        if "DEBUG_REG" not in os.environ and not user_config["model_dir"]:
            if user_config["repo_type"] == "git":
                fp.write("Commits from last day:\n")
                with open(results_dir + "/gitLog", "r") as inf:
                    fp.write(inf.read())
            fp.write("\n")
            fp.write("-" * 80 + "\n")
        if "yes" in use_html:
            fp.write("</pre><html>\n")

    subject = '"[NU-WRF-regression]" '
    cmd = "/usr/bin/mail -s " + subject + mail_to + " < " + test_report
    if "yes" in use_html:
        pref = 'mutt -e "set content_type=text/html" -s '
        cmd = pref + subject + mail_to + " < " + test_report
    sp.call(cmd, shell=True)


def get_compilers(config):
    comp_config = utils.config_section_map(config, "COMPCONFIG")
    comp_list = comp_config["compilers"].split(",")
    compilers = []
    [compilers.append(c.strip()) for c in comp_list]
    return compilers


def get_modulelist(config):
    comp_config = utils.config_section_map(config, "COMPCONFIG")
    m_list = comp_config["modulelist"].split(",")
    module_list = []
    [module_list.append(m.strip()) for m in m_list]
    return module_list
