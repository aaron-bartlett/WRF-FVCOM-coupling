"""
  This script executes NU-WRF run combinations and generates
  results used to verify various reproducibility measures.
"""
import shared.utils as utils
import logging
import RegRuns
import sys


def main():
    test_case = sys.argv[1]

    logger = logging.getLogger("MAIN    ")

    run = RegRuns.RegRuns(test_case)
    # Setup a logging stream for this run
    run.set_logging()

    report_file = (
        run.get_opt("scratch_dir")
        + "/results/"
        + "/"
        + run.get_opt("compiler")
        + "."
        + run.get_opt("time_stamp")
        + "/"
        + run.name
        + ".diff"
    )
    with open(report_file, "w") as fp:
        logger.info(
            "<"
            + run.name
            + "> verification type: <"
            + run.verification
            + ">, compiler: <"
            + run.compiler
            + ">"
        )
        for component in run.components.split(","):
            rc = run.run(component)
            # rc is the return code from "grep", rc=0 means failure
            if rc == 0:
                break
        if rc == 1:
            logger.info("Checking run output against baseline")
            RegRuns.check_baseline(run)
        utils.write_results(run.results, fp)

    logger.info(run.name + " verification complete.")


if __name__ == "__main__":
    main()
