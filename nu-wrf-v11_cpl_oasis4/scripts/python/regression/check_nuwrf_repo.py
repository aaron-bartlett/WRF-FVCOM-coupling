#!/usr/local/other/python/GEOSpyD/2019.10_py3.7/2020-01-15/bin/python
import os
import subprocess


def syscmd(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    p = subprocess.Popen(cmd, shell=True, stdout=stdout, stderr=stderr)
    (out, err) = p.communicate()
    return out.decode('utf-8')


def main():
    print("Start...")
    #branches = {"master": "NO"}
    branches = {"nompi": "NO"}
    #root = "/discover/nobackup/ccruz/devel/nu-wrf/code/gitlab/nu-wrf-dev"
    root = "/discover/nobackup/ccruz/devel/nu-wrf/code/gitlab/"
    python_cmd = "/usr/local/other/python/GEOSpyD/2019.10_py3.7/2020-01-15/bin/python"
    git_exe = "/usr/local/other/git/2.24.0/libexec/git-core/git"
    old_path = os.environ["PATH"]
    os.environ["PATH"] = git_exe + os.pathsep + old_path

    print("Check for repository updates...")
    for branch, dotest in branches.items():
        os.chdir(root + branch)
        before = syscmd(git_exe+" rev-parse HEAD")
        rc = syscmd(git_exe+" pull")  # update repository for "next" time
        after = syscmd(git_exe+" rev-parse HEAD")
        if str(before) != str(after):
            print(" -- Changes detected in branch " + branch)
            branches[branch] = "YES"
        else:
            print(" -- NO changes detected in branch " + branch)

    for branch, dotest in branches.items():
        if "YES" in dotest:
            print(" -- Running tests for branch " + branch)
            os.chdir(root + branch + "/scripts/python/regression/")
            with open(branch + ".cron", "w") as log:
                out = syscmd("module avail", stdout=log, stderr=log)
                #out = syscmd(python_cmd + " reg " + branch, stdout=log, stderr=log)
        else:
            with open(branch + ".cron", "w") as log:
                subprocess.call("module avail", shell=True)
            #print(" -- Nothing to do for branch " + branch)


if __name__ == "__main__":
    main()
