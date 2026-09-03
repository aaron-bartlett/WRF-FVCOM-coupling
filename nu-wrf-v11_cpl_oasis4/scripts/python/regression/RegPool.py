import time
import subprocess
import shlex
import multiprocessing
import logging

logger = logging.getLogger("pool")

# Reference: https://pymotw.com/2/multiprocessing/communication.html


class Worker(multiprocessing.Process):
    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                logger.debug("%s: Exiting" % proc_name)
                self.task_queue.task_done()
                break
            logger.debug("%s: %s" % (proc_name, next_task))
            answer = next_task()
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return


class Task(object):
    def __init__(self, a):
        self.a = a

    def __call__(self):
        process = subprocess.Popen(self.a, shell=True)
        while process.poll() is None:
            time.sleep(10)
        if process.returncode != 0:
            logger.error("%r failed: %s" % (self.a, process.returncode))
        logger.debug("%r is done" % self.a)

    def __str__(self):
        return "%s" % self.a


class Batch(object):
    jobs = []

    def __init__(self, cmd):
        self.cmd = cmd

    def __call__(self):
        self.slurm_command()
        while True:
            # If list is empty then we are done
            if len(self.jobs) == 0:
                break
            for job in self.jobs:
                rc = subproc("squeue -j " + str(job) + " -t PD,R -h -o %t")
                # If job is done, remove from list
                if not rc:
                    logger.debug("..." + job + " is done")
                    self.jobs.remove(job)
                else:
                    if "PD" in rc:
                        logger.debug("..." + job + " is pending")
                    elif "R" in rc:
                        logger.debug("..." + job + " is running")
                    elif "U" in rc:
                        logger.debug("..." + job + " is terminating")
                    else:
                        pass
                    time.sleep(30)

    def __str__(self):
        return "%s" % self.cmd

    # This function submits a batch job under SLURM and creates a jobs list
    # NCCS-DISCOVER only.
    def slurm_command(self):
        rc = subproc(self.cmd % vars())
        logger.debug(rc)
        # output example: "Submitted batch job 12345", so we
        # parse the output from rc and append job ID to jobs list
        self.jobs.append(shlex.split(rc)[3])


def subproc(cmd):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    return out.decode("utf-8")


def run_commands(cmds, use_batch):
    # Establish communication queues
    commands = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    # Start workers
    num_jobs = len(cmds)
    logger.debug("Creating %d workers" % num_jobs)
    workers = [Worker(commands, results) for _ in range(num_jobs)]
    [w.start() for w in workers]

    # Enqueue jobs
    for cmd in cmds:
        if use_batch == "yes":
            commands.put(Batch(cmd))
        else:
            commands.put(Task(cmd))

    # add one stop value per worker to the job queue
    [commands.put(None) for _ in range(num_jobs)]

    # Wait for all of the commands to finish
    commands.join()
