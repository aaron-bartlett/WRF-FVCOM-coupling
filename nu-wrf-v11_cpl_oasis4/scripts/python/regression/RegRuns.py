"""
Class to instantiate NU-WRF regression testing runs.
"""
from __future__ import print_function
import sys
import glob
import os
import subprocess as sp
import logging

try:
    import ConfigParser
except ImportError:
    from configparser import ConfigParser
import shlex
import shutil
import shared.utils as utils
import RegTest


class RegRuns(RegTest.RegTest):
    def __init__(self, run_name):
        super(RegRuns, self).__init__(run_name)

        # Override defaults with values specified in config file
        cfg = get_configuration(self)
        self.set_opts(cfg, self.name)
        # In case we want to debug this run:
        if "DEBUG_REG" in os.environ:
            self.debug = True

    # Setup a logging object for each run. Note that the output file gets
    # logging DEBUG while STDOUT only gets logging INFO in order to minimize
    # verbosity.
    def set_logging(self):
        # Note filemode is 'append' because we want one logger per job
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            datefmt="%m-%d %H:%M",
            filename=self.get_opt("scratch_dir")
            + "/results/"
            + "/"
            + self.compiler
            + "."
            + self.time_stamp
            + "/"
            + self.name
            + "-regression.log",
            filemode="a",
        )
        stdout_log = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter("%(name)s : %(message)s")
        stdout_log.setFormatter(formatter)
        stdout_log.setLevel(logging.INFO)
        logger = logging.getLogger()
        logger.addHandler(stdout_log)

    def build(self, command):
        logger = logging.getLogger("BLDCOMP ")
        logger.info(command + " build: <" + self.compiler + ">")
        cmd = "./" + command
        return self.run_shell_command(cmd)

    def run(self, component):
        logger = logging.getLogger("RUNCOMP ")
        logger.info(component + " run: <" + self.compiler + ">")
        cmd = "./" + component + ".reg"
        return self.run_shell_command(cmd)

    def run_shell_command(self, command_line):
        # Function used to run scripts
        command_line_args = shlex.split(command_line)
        logger = logging.getLogger("SYSTEM  ")
        logger.info('Subprocess: "' + command_line + '"')
        if self.debug:
            logger.info(command_line)
            self.results = [
                self.name,
                self.compiler.center(10, " "),
                "r+".center(10, " "),
                "-".center(10, " "),
            ]
            return 0
        else:
            ok = 1 # success
            try:
                command_line_process = sp.Popen(
                    command_line_args, stdout=sp.PIPE, stderr=sp.STDOUT
                )
                command_line_process.communicate()
            # Very unlikely this will happen:
            except OSError as e:
                sys.exit("failed to execute command: %s" % (str(e)))

            # Parse file output from components
            finally:
                logfile = None
                if "./geogrid" in command_line:
                    logfile = "geogrid_logs/geogrid*"
                    grepstr = "Successful"
                elif "./ungrib" in command_line:
                    logfile = "ungrib_logs/ungrib.log"
                    grepstr = "Successful"
                elif "./metgrid" in command_line:
                    logfile = "metgrid_logs/metgrid*"
                    grepstr = "Successful"
                elif "./real" in command_line:
                    logfile = "real_logs/real.rsl.out*"
                    grepstr = "SUCCESS"
                elif "./ideal" in command_line:
                    logfile = "ideal_logs/rsl.out*"
                    grepstr = "SUCCESS"
                elif "./wrf" in command_line:
                    logfile = "wrf_logs/wrf.rsl.out*"
                    grepstr = "SUCCESS"
                elif "./lis." in command_line:
                    logfile = "lis_logs/lislog*"
                    grepstr = "LIS Run completed"
                elif "./ldt_prelis" in command_line:
                    logfile = "ldtlog_prelis.0000"
                    grepstr = "Finished LDT run"
                elif "./ldt_postlis" in command_line:
                    logfile = "ldtlog_postlis.0000"
                    grepstr = "Finished LDT run"
                else:
                    logfile = ""

                if logfile:
                    logfiles = glob.glob(logfile)
                    # An ugly hack needed for "special" cases:
                    if "ldt_prelis" in command_line:
                        logfiles = ["ldtlog_prelis.0000"]
                    elif "ldt_postlis" in command_line:
                        logfiles = ["ldtlog_postlis.0000"]
                    if logfiles:
                        grep = sp.Popen(
                            ["grep", grepstr, logfiles[0]], stdout=sp.PIPE, stderr=sp.PIPE
                        )
                        out, err = grep.communicate()
                        print(out.decode("utf-8"))
                        if out:
                            logger.info(out.decode("utf-8"))
                    else:
                        ok = 0

                if ok == 1:
                    self.results = [
                        self.name,
                        self.compiler.center(10, " "),
                        "r+".center(10, " "),
                        "-".center(10, " "),
                    ]
                    logger.info(command_line + " SUCCESS")
                else:
                    logger.error(command_line + " FAILURE")
                    self.results = [
                        self.name,
                        self.compiler.center(10, " "),
                        "rF".center(10, " "),
                        "-".center(10, " "),
                    ]
            return ok

def get_configuration(run: object) -> object:
    configfile = "./" + run.name + ".cfg"
    if os.path.isfile(configfile):
        try:
            config = ConfigParser.ConfigParser()
        except AttributeError:
            config = ConfigParser()
        config.read(configfile)
        run.parser_to_test(config)
    else:
        print(" *** Configuration file does not exist *** " + configfile)
        sys.exit(1)
    return config


def check_baseline(run):
    logger = logging.getLogger("BASE_VER")
    base_dir = run.get_opt("baseline_dir")
    compiler = run.get_opt("compiler")

    diff = run.get_opt("scripts_dir") + "/diffnc4.py "
    dest_dir = base_dir + "/" + compiler + "/" + run.name

    logger.debug("Begin verification.")

    if "casa" in run.name:
        outfiles = glob.glob("CO2_d01*")
    elif "nldas2" in run.name:
        outfiles = glob.glob("LIS_HIST*")
    else:
        outfiles = glob.glob("wrfout*")
    # If no outfiles, run must have failed.
    if len(outfiles) == 0:
        logger.debug("No output files, returning...")
        run.results[2] = "rF"
        return

    # If expected_output is not specified then simply update files
    if not run.get_opt("expected_output"):
        logger.debug("expected_output was not specified.")
        # Make sure we really want to update destination
        if run.update_base == "yes":
            utils.mkdir_p(dest_dir)
            for f in outfiles:
                logger.debug("Update baseline " + f)
                shutil.copy(f, dest_dir + "/" + f)
                run.results[3] = "C"
        else:
            run.results[3] = "vF"
        return

    # expected_output is specified. Check for ALL expected output
    for f in run.get_opt("expected_output").split(","):
        # If ANY expected output file is missing, run must have failed
        # or perhaps the checkpoint frequency was changed (user mistake).
        if not os.path.isfile(f):
            logger.error("The specified expected file [" + f + "] was not found!")
            run.results[2] = "rF"
            return

    # Ok, ALL the expected output is available. Let's compare with baseline
    utils.mkdir_p(dest_dir)
    for f in run.get_opt("expected_output").split(","):
        dest_file = dest_dir + "/" + f
        cmd = diff + f + " " + dest_file
        logger.debug(cmd)
        out = check_cmp(cmd)
        # Comparison (diff) yielded nothing. That means zero-diff:
        if not out:
            logger.debug("Zero-diff")
            run.results[3] = "+"
        else:
            if not os.path.isfile(dest_file):
                logger.debug("Creating baseline file " + f)
                shutil.copy(f, dest_file)
                run.results[3] = "C"
            else:
                if run.update_base == "yes":
                    logger.debug("Update baseline " + f)
                    shutil.copy(f, dest_file)
                run.results[3] = "vF"


def check_cmp(command_line):
    command_line_args = shlex.split(command_line)
    command_line_process = sp.Popen(command_line_args, stdout=sp.PIPE, stderr=sp.STDOUT)
    command_line_process.wait()
    stdout, stderr = command_line_process.communicate()
    return stdout.decode("utf-8")
