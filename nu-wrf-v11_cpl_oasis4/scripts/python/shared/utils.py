"""
General utilities for NU-WRF building and testing
"""
import os
import sys
import errno
import datetime
import shutil
import shlex
import re
import tempfile
import logging
import io

try:
    import ConfigParser
except ImportError:
    import configparser
import subprocess as sp

logger = logging.getLogger("utils")


def replace_infile(filename, src, dst):
    with open(filename, "r") as f:
        lines = f.readlines()
        with open(filename, "w") as f:
            for line in lines:
                if src in line:
                    line = line.replace(src, dst)
                f.write(line)


def textio_seek(fobj, amount, whence=0):
    fobj.buffer.seek(amount, whence)
    return io.TextIOWrapper(
        fobj.buffer, encoding=fobj.encoding, errors=fobj.errors, newline=fobj.newlines
    )


def set_time_stamp():
    return datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def get_config(filename):
    if not os.path.isfile(filename):
        logger.error(" +++ Error: " + filename + " : file does not exist +++")
        sys.exit(1)
    try:
        config = ConfigParser.SafeConfigParser()
    except NameError:
        config = configparser.ConfigParser(
            interpolation=configparser.ExtendedInterpolation()
        )
    config.read(filename)
    
    logger.debug("Done reading %s", filename)
    return config


def config_section_map(config, section):
    # Return a dict (i.e. a key,value pair) from each section in config file
    section_dict = {}
    options = config.options(section)

    for option in options:
        try:
            section_dict[option] = config.get(section, option)
            if section_dict[option] == -1:
                logger.info("skip option %s", option)
        except KeyError:
            logger.error("exception on %s!" % option)
            section_dict[option] = None

    return section_dict


def clean_scratch(config):
    logger.debug("Clean up testing environment...")
    user_config = config_section_map(config, "USERCONFIG")
    results_dir = user_config["scratch_dir"] + "/results"
    builds_dir = user_config["scratch_dir"] + "/builds"

    if not os.path.exists(results_dir):
        mkdir_p(results_dir)
    else:
        clean_dir(results_dir)

    if not os.path.exists(builds_dir):
        mkdir_p(builds_dir)
    else:
        clean_dir(builds_dir)


def get_repotype(url):
    if url[0:3] == "svn" or url[0:4] == "file" or url[0:4] == "http":
        return "svn"
    else:
        return "git"

def force_symlink(src, dst):
    try:
        os.symlink(src, dst)
    except FileExistsError as e:
        if e.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src, dst)

def mkdir_p(path):
    # Layer on top of os.makedirs - with error checking
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            logger.error("Permission denied: cannot create " + path)
            sys.exit()


def clean_dir(adir):
    # "Safe" way to clean the contents of a directory
    if adir == "/" or adir == "\\":
        logger.error("Cannot clean %s", adir)
        return
    else:
        for file_object in os.listdir(adir):
            logger.debug("Will clean up %s", adir)
            file_object_path = os.path.join(adir, file_object)
            if os.path.isfile(file_object_path):
                os.unlink(file_object_path)
            else:
                try:
                    shutil.rmtree(file_object_path)
                except OSError:
                    logger.error("Permission denied: cannot remove " + file_object_path)
                    sys.exit()


def copy_file(src, dest):
    try:
        shutil.copy(src, dest)
        return 0
    # eg. src and dest are the same file
    except shutil.Error as e:
        logger.error("Error: %s" % e)
    # eg. source or destination doesn't exist
    except IOError as e:
        logger.error("Error: %s" % e)
    return 1


def which(program):
    # Test if an executable program exists in the path - like unix's which
    def is_exe(file_path):
        return os.path.isfile(file_path) and os.access(file_path, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def logger_setup(logger_name):
    # log messages to a file
    logging.basicConfig(
        filename=str(logger_name) + ".LOG",
        format="%(levelname)s %(module)s:%(lineno)s %(funcName)s %(message)s",
        level=logging.DEBUG,
        filemode="w",
    )
    # Add a stream handler (using stdout) to default logger
    stdout_log = logging.StreamHandler(sys.stdout)
    # Set level to INFO
    stdout_log.setLevel(logging.INFO)
    formatter = logging.Formatter("%(levelname)s %(name)s : %(message)s")
    stdout_log.setFormatter(formatter)
    root = logging.getLogger()
    root.addHandler(stdout_log)


def show_config(filename):
    try:
        config = ConfigParser.SafeConfigParser()
    except NameError:
        config = configparser.ConfigParser()
    config.read(filename)

    output = sp.run(["git", "rev-parse", "--show-toplevel"], stdout=sp.PIPE)
    repo = output.stdout.decode("utf-8").rstrip()
    for section_name in config.sections():
        match = not re.search("CONFIG", section_name)
        if match:
            print("-" * 80)
            print("TESTCASE:", section_name)
            print("-" * 80)
            tokens = section_name.split("_", 1)
            readme = repo + "/testcases/" + tokens[0] + "/" + tokens[1] + "/README"
            try:
                with open(readme, "r") as f:
                    print(f.read())
            except FileNotFoundError:
                print()
            print("Settings:")
            for name, value in config.items(section_name):
                print("  %s = %s" % (name, value))
            print()


def write_results(results, fp):
    fp.write("{0:^40}".format(results[0]))
    fp.write("{0:^20}".format(results[1]))
    fp.write("{0:^10}".format(results[2]))
    fp.write("{0:^10}".format(results[3]))
    fp.write("\n")


def header(cfg):
    logger.info("USER CONFIGURATION:")
    user_config = config_section_map(cfg, "USERCONFIG")
    mail_to = user_config["mail_to"]
    branch = user_config["repo_branch"]
    if not branch:
        branch = "none specified"
    repo_type = user_config["repo_type"]
    build_type = user_config["build_type"]
    sections = cfg.sections()

    print("-" * 80)
    if user_config["exe_dir"]:
        print("Using executables in: " + user_config["exe_dir"])
    else:
        print("Testing _" + repo_type + "_ repository: " + user_config["repo_url"])
        print("Build type: _" + build_type + "_, branch: _" + branch + "_")
    print("Output in: " + user_config["scratch_dir"])
    print("Configurations (verification type):")
    for section_name in sections:
        if "CONFIG" not in section_name:
            temp_dict = dict(cfg.items(section_name))
            printed = False  # avoid printing section_name more than once
            for name, value in list(temp_dict.items()):
                if "verification" in name and not printed:
                    print("\t", section_name, " (", value, ")")
                    printed = True
                else:
                    if not printed:
                        print("\t", section_name, " (run)")
                        printed = True
    if mail_to:
        print("Results will be mailed to: " + mail_to)
    print("-" * 80)


def grep_file(file, pattern):
    with open(file, "r") as f:
        for line in f:
            if pattern in line:
                return True
        return False


def find_files(in_dir, pattern):
    # Use like Unix's find
    import fnmatch

    file_list = []
    # Walk through directory
    for dName, sdName, fList in os.walk(in_dir):
        for fileName in fList:
            if fnmatch.fnmatch(fileName, pattern):  # Match search string
                file_list.append(os.path.join(dName, fileName))
    return file_list


def sed_inplace(filename, pattern, repl):
    # https://stackoverflow.com/questions/4427542/how-to-do-sed-like-text-replace-with-python
    # Perform the pure-Python equivalent of in-place `sed` substitution: e.g.,
    # `sed -i -e 's/'${pattern}'/'${repl}' "${filename}"`.
    # For efficiency, precompile the passed regular expression.
    pattern_compiled = re.compile(pattern)

    # For portability, NamedTemporaryFile() defaults to mode "w+b" (i.e., binary
    # writing with updating). This is usually a good thing. In this case,
    # however, binary writing imposes non-trivial encoding constraints trivially
    # resolved by switching to text writing. Let's do that.
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
        with open(filename) as src_file:
            for line in src_file:
                tmp_file.write(pattern_compiled.sub(repl, line))

    # Overwrite the original file with the munged temporary file in a
    # manner preserving file attributes (e.g., permissions).
    shutil.copystat(filename, tmp_file.name)
    shutil.move(tmp_file.name, filename)


def sp_call(cmd):
    # Run an OS command
    if "DEBUG" in os.environ:
        logger.debug("--- COMMAND: " + cmd)
        return
    try:
        sp.call(cmd, shell=True)
    except sp.CalledProcessError as ex:
        logger.error("Error {} in command".format(ex.returncode))


def sp_check_call_make_log(command_line, envs=None):
    command_line_args = shlex.split(command_line)
    # If used, this function will capture stdout/stderr in one make.log file
    if "DEBUG" in os.environ:
        logger.debug("--- COMMAND: " + command_line)
        return
    try:
        fout = open(os.environ.get("NUWRFDIR") + "/make.log", "a")
        sp.check_call(command_line_args, env=envs, stdout=fout, stderr=sp.STDOUT)
        fout.close()
    # Very unlikely this will happen:
    except OSError as e:
        sys.exit("failed to run shell: %s" % (str(e)))
    return 0


def run_shell_command(command_line, envs=None):
    command_line_args = shlex.split(command_line)
    logger.debug('Subprocess: "' + command_line + '"')
    if "DEBUG" in os.environ:
        logger.debug("--- COMMAND: " + command_line)
        return

    else:
        try:
            command_line_process = sp.Popen(
                command_line_args, env=envs, stdout=sp.PIPE, stderr=sp.PIPE
            )
            stdout, stderr = command_line_process.communicate()
            if os.getenv("DEBUG_NUWRF_SHELL_CMD"):
                logger.error(stderr.decode("utf-8"))
        # Very unlikely this will happen:
        except OSError as e:
            logger.error("failed to execute command: %s" % (str(e)))
            sys.exit()
        except sp.CalledProcessError as err:
            logger.error("Error {} in command".format(err))
        else:
            # Set DEBUG_WRF to print STDOUT
            if "DEBUG_NUWRF" in os.environ:
                logger.info("--- COMMAND: " + command_line)
                logger.info(stdout.decode("utf-8"))
        # Check return code from script
        finally:
            return command_line_process.returncode


def run_configure(command_line, heredoc=None, envs=None):
    command_line_args = shlex.split(command_line)
    logger.debug('Subprocess: "' + command_line + '"')
    if "DEBUG" in os.environ:
        logger.debug("--- COMMAND: " + command_line)
        return

    else:
        try:
            command_line_process = sp.Popen(
                command_line_args,
                env=envs,
                universal_newlines=True,
                stdin=sp.PIPE,
                stdout=sp.PIPE,
                stderr=sp.PIPE
            )
            command_line_process.stdin.write(heredoc)
            command_line_process.stdin.flush()
            command_line_process.wait()
        # Very unlikely this will happen:
        except OSError as e:
            logger.error("failed to execute command: %s" % (str(e)))
            sys.exit()
        except sp.CalledProcessError as err:
            logger.error("Error {} in command".format(err))
        finally:
            if "DEBUG_NUWRF" in os.environ:
                for line in iter(command_line_process.stdout.readline, ''):
                    print(line); sys.stdout.flush()
            return command_line_process.returncode

