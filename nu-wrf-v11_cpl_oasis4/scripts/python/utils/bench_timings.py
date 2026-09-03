#!/usr/bin/env python
# ------------------------------------------------------------------------------
# NASA/GSFC, CISTO, Code 606, Advanced Software Technology Group
# ------------------------------------------------------------------------------
#
# MODULE: bench_timings
#
# AUTHOR:
# Eric Kemp, NASA CISTO/SSAI
#
# DESCRIPTION:
# Parses RSL standard error file from WRF model to calculate breakdown of
# time spent in different WRF components called by the solver. WRF must have
# been compiled with the (unadvertised) -DBENCH flag passed to the cpp
# preprocessor program.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Import Python standard modules
# ------------------------------------------------------------------------------

from __future__ import print_function
import operator
import sys

# ------------------------------------------------------------------------------
# Get name of RSL standard error file to process.
# ------------------------------------------------------------------------------

if len(sys.argv) != 2:
    text = "ERROR, invalid command line arguments!\n"
    text += "Usage: %s rslErrorFileName\n" % (sys.argv[0])
    text += "  where rslErrorFileName is RSL standard error file from WRF,\n"
    text += "  e.g., rsl.error.0000"
    raise SystemExit(text)
rslfile = sys.argv[1]

# ------------------------------------------------------------------------------
# Parse the RSL standard error file, summing the microseconds recorded in each
# WRF component in each computational time step. The sums are stored in a
# Python dictionary with keys corresponding to each WRF component where data
# are provided.
#
# The RSL file should include output like this:
#
# Timing for main: time 2009-04-10_12:00:36 on domain   1:   57.29905 elapsed seconds
# solve_tim=       115396
# step_prep_tim=         1044
# set_phys_bc_tim=          384
# ...etc...
# ------------------------------------------------------------------------------

timings = {}
lines = open(rslfile, "r").readlines()
for line in lines:
    firstWord = line.split()[0]
    if "_tim=" in firstWord:
        if firstWord not in list(timings.keys()):
            timings[firstWord] = 0
        microseconds = int(line.split()[1])
        timings[firstWord] += microseconds

# ------------------------------------------------------------------------------
# "solve_tim" is the total time spent in the solver. All the other components
# are subsets.
# ------------------------------------------------------------------------------

total = timings["solve_tim="]
print("Total time: ", total, " microseconds")

# ------------------------------------------------------------------------------
# Now output each component relative and absolute time. But sort the components
# in descending order.
# ------------------------------------------------------------------------------

sorted_timings = sorted(list(timings.items()), key=operator.itemgetter(1))
sorted_timings.reverse()
sum = 0
for (key, microseconds) in sorted_timings:
    if "solve_tim" in key:
        print(
            "%s : ( %6.2f%% ) ( %8d microseconds )"
            % (key[:-1], float(microseconds) / float(total) * 100.0, microseconds)
        )
    else:
        print(
            "\t%23s : ( %6.2f%% ) ( %8d microseconds )"
            % (key[:-1], float(microseconds) / float(total) * 100.0, microseconds)
        )
        sum += microseconds
