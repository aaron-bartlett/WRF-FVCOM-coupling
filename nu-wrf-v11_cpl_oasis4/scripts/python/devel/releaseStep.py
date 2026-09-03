#!/usr/bin/env python
import subprocess
import sys
import os
import shutil


def execute(cmd):
    try:
        p = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        if e.returncode == 127:
            sys.exit("returned code {}".format(e.returncode))
        elif e.returncode <= 125:
            sys.exit("command failed, returned code {}".format(e.returncode))
        else:
            sys.exit("command crashed, shell returned code {}".format(e.returncode))
    except OSError as e:
        sys.exit("failed to run shell: {}".format(str(e)))


if len(sys.argv) < 2:
    print("Usage: " + sys.argv[0] + " version_string (e.g. v10-wrf391-lisf)")
    sys.exit()
version = sys.argv[1]

top = os.getcwd()
url = "git@git.mysmce.com:astg/nuwrf/nu-wrf-dev.git"
prefix = "nu-wrf_"
dirs = [
    "/discover/nobackup/projects/nu-wrf/releases/stable",
    "/nobackupp8/nuwrf/releases/stable",
]
tarballs = [prefix + version + ".bz2", prefix + version + ".tgz"]
branch = "develop"

print("Cloning " + url + "...")
execute("git clone -b " + branch + " " + url + " nu-wrf-dev")
os.chdir("nu-wrf-dev")
for f in tarballs:
    print("Archiving into " + f + "...")
    compress = "gzip" if "tgz" in f else "bzip2"
    cmd = (
        "git archive --format=tar --prefix=nu-wrf-dev/ " + branch + "|" + compress + ">" + f
    )
    execute(cmd)
    execute("chgrp s0942 " + f)
    execute("chmod 640 " + f)
    print("Releasing...")
    print("...discover")
    execute("cp -p " + f + " " + dirs[0])
    #print("...pleiades")
    #execute("rsync -av " + f + " pfe.nas.nasa.gov:" + dirs[1])

print("Cleaning...")
os.chdir(top)
shutil.rmtree(top + "/nu-wrf-dev")
