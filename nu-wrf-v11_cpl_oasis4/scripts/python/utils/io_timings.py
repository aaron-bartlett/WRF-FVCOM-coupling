#!/usr/bin/env python
# ------------------------------------------------------------------------------
# NASA/GSFC, CISTO, Code 606, Advanced Software Technology Group
# ------------------------------------------------------------------------------
#
# MODULE: io_timings
#
# AUTHOR:
# Eric Kemp, NASA CISTO/SSAI
#
# DESCRIPTION:
# Parses RSL output file from WRF model to calculate breakdown of
# time spent in computations and I/O.
#
# This script does not consider the breakdown of runtime in each part of
# the WRF solver. That is handled by bench_timings.py
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
# Parse the RSL output file. summing the elapsed seconds recorded for each
# time step and I/O point. The sums are stored in a Python dictionary with
# keys corresponding to each process where data are provided.
#
# The RSL files should include output like this:
#
# Timing for main: time 2009-04-10_14:58:48 on domain   1:    0.10355 elapsed seconds
# Timing for main: time 2009-04-10_14:59:24 on domain   1:    0.10242 elapsed seconds
# Timing for main: time 2009-04-10_15:00:00 on domain   1:    0.10422 elapsed seconds
# Timing for Writing wrfout_d01_2009-04-10_15:00:00 for domain        1:    4.18171 elapsed seconds
# Timing for Writing restart for domain        1:   13.55360 elapsed seconds
# Timing for processing lateral boundary for domain        1:    0.20064 elapsed seconds
# Timing for main: time 2009-04-10_15:00:36 on domain   1:   18.35575 elapsed seconds
# ...etc...
# ------------------------------------------------------------------------------

# Individual restart file names are not in RSL file, so we have to do some
# statistics here.
writingAllRestartsMax = 0
writingAllRestartsMin = 99999999999999
writingAllRestartsMean = 0
writingAllRestartsCount = 0

timings = {}
lines = open(rslfile, "r").readlines()
for line in lines:
    if "elapsed seconds" in line:
        word = line.split()[2]
        if word == "main:":
            key = "main:"
        else:
            key = "%s %s" % (word, line.split()[3])
        if key not in list(timings.keys()):
            timings[key] = 0
        seconds = float(line.split()[-3])
        timings[key] += seconds
        if key == "Writing restart":
            writingAllRestartsMin = min(writingAllRestartsMin, seconds)
            writingAllRestartsMax = max(writingAllRestartsMax, seconds)
            writingAllRestartsCount += 1

writingAllRestartsTotal = timings["Writing restart"]
timings["AllWrfRestartsWriting"] = writingAllRestartsTotal
if writingAllRestartsCount > 0:
    writingAllRestartsMean = writingAllRestartsTotal / float(writingAllRestartsCount)

# ------------------------------------------------------------------------------
# Add up all the writing times.
# ------------------------------------------------------------------------------

writingAllMax = 0
writingAllMin = 99999999999999
writingAllMean = 0
total = 0
count = 0
for key in list(timings.keys()):
    if "Writing" in key:
        total += timings[key]
        count += 1
        writingAllMin = min(writingAllMin, timings[key])
        writingAllMax = max(writingAllMax, timings[key])
timings["AllWriting"] = total
if count > 0:
    writingAllMean = total / float(count)
writingAllCount = count

writingAllPressMax = 0
writingAllPressMin = 99999999999999
writingAllPressMean = 0
total = 0
count = 0
for key in list(timings.keys()):
    if "Writing wrfpress" in key:
        total += timings[key]
        count += 1
        writingAllPressMin = min(writingAllPressMin, timings[key])
        writingAllPressMax = max(writingAllPressMax, timings[key])
timings["AllWrfPressWriting"] = total
if count > 0:
    writingAllPressMean = total / float(count)
writingAllPressCount = count

writingAllOutMax = 0
writingAllOutMin = 99999999999999
writingAllOutMean = 0
total = 0
count = 0
for key in list(timings.keys()):
    if "Writing wrfout" in key:
        total += timings[key]
        count += 1
        writingAllOutMin = min(writingAllOutMin, timings[key])
        writingAllOutMax = max(writingAllOutMax, timings[key])
timings["AllWrfOutWriting"] = total
if count > 0:
    writingAllOutMean = total / float(count)
writingAllOutCount = count

writingAll2dOutMax = 0
writingAll2dOutMin = 99999999999999
writingAll2dOutMean = 0
total = 0
count = 0
for key in list(timings.keys()):
    if "Writing wrf2dout" in key:
        total += timings[key]
        count += 1
        writingAll2dOutMin = min(writingAll2dOutMin, timings[key])
        writingAll2dOutMax = max(writingAll2dOutMax, timings[key])
timings["AllWrf2dOutWriting"] = total
if count > 0:
    writingAll2dOutMean = total / float(count)
writingAll2dOutCount = count

writingAllDiagnosticsMax = 0
writingAllDiagnosticsMin = 99999999999999
writingAllDiagnosticsMean = 0
total = 0
count = 0
for key in list(timings.keys()):
    if "Writing wrfdiagnostics" in key:
        total += timings[key]
        count += 1
        writingAllDiagnosticsMin = min(writingAllDiagnosticsMin, timings[key])
        writingAllDiagnosticsMax = max(writingAllDiagnosticsMax, timings[key])
timings["AllWrfDiagnosticsWriting"] = total
if count > 0:
    writingAllDiagnosticsMean = total / float(count)
writingAllDiagnosticsCount = count

# ------------------------------------------------------------------------------
# "main" is the total time spent in each time step. All other components are
# subsets.
# ------------------------------------------------------------------------------

total = timings["main:"]
print("Total time: ", total, " seconds")

# ------------------------------------------------------------------------------
# Now output each I/O component relative and absolute time. But sort the
# components in descending order. Also provide averages for "All Write"
# cases.
# ------------------------------------------------------------------------------

sorted_timings = sorted(list(timings.items()), key=operator.itemgetter(1))
sorted_timings.reverse()
for (key, seconds) in sorted_timings:
    if key not in [
        "main:",
        "AllWriting",
        "AllWrfPressWriting",
        "AllWrfOutWriting",
        "AllWrf2dOutWriting",
        "AllWrfDiagnosticsWriting",
        "AllWrfRestartsWriting",
    ]:
        continue

    print("%s" % (key))
    print(
        "( Pct of run: %6.2f%% ) ( %.2f sec )"
        % (float(seconds) / float(total) * 100.0, seconds)
    )

    if key == "AllWriting":
        print("( number of writes: %d )" % (writingAllCount))
        print(
            "( mean/max/min sec per write: %.2f %.2f %.2f )"
            % (writingAllMean, writingAllMax, writingAllMin)
        )
    elif key == "AllWrfPressWriting":
        print("( number of writes: %d )" % (writingAllPressCount))
        print(
            "( mean/max/min sec per write: %.2f %.2f %.2f )"
            % (writingAllPressMean, writingAllPressMax, writingAllPressMin)
        )
    elif key == "AllWrfOutWriting":
        print("( number of writes: %d )" % (writingAllOutCount))
        print(
            "( mean/max/min sec per write: %.2f %.2f %.2f )"
            % (writingAllOutMean, writingAllOutMax, writingAllOutMin)
        )
    elif key == "AllWrf2dOutWriting":
        print("( number of writes: %d )" % (writingAll2dOutCount))
        print(
            "( mean/max/min sec per write: %.2f %.2f %.2f )"
            % (writingAll2dOutMean, writingAll2dOutMax, writingAll2dOutMin)
        )
    elif key == "AllWrfDiagnosticsWriting":
        print("( number of writes: %d )" % (writingAllDiagnosticsCount))
        print(
            "( mean/max/min sec per write: %.2f %.2f %.2f)"
            % (
                writingAllDiagnosticsMean,
                writingAllDiagnosticsMax,
                writingAllDiagnosticsMin,
            )
        )
    elif key == "AllWrfRestartsWriting":
        print("( number of writes: %d )" % (writingAllRestartsCount))
        print(
            "( mean/max/min sec per write: %.2f %.2f %.2f)"
            % (writingAllRestartsMean, writingAllRestartsMax, writingAllRestartsMin)
        )
