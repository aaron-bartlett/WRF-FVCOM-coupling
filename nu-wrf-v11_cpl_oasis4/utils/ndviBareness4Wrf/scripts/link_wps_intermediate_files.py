#!/usr/bin/env python
#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Sciences and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# SCRIPT: link_wps_intermediate_files.py
#
# AUTHOR:
# Eric Kemp, NASA GSFC/SSAI
#
# DESCRIPTION:
# Creates symbolic links to WPS intermediate files.  Although the WPS format 
# includes internal time records for individual fields, METGRID does not use 
# them, and instead assumes the valid date/time matches the date/time used in 
# the file name.  Thus symbolic links can be used to force METGRID to include 
# these fields in time-varying lateral boundary conditions without wasting 
# disk space.
#
# REVISION HISTORY:
# 01 Sep 2016 - Initial version.
#
#------------------------------------------------------------------------------

# Standard Python modules
import ConfigParser
import datetime
import glob
import optparse
import os
import sys

# Useful constants
SECTION_NAME = "LINK_WPS_INT"

# Check command line arguments
parser = optparse.OptionParser(usage="%prog [--help] [-P]")
parser.add_option("-P","--prefix",dest="prefix",
                  default=None,
                  help="Override prefix setting in config file")
(options,args) = parser.parse_args()

# Read config file
config = ConfigParser.RawConfigParser()
config.read("link_wps_int.cfg")
prefix = config.get(SECTION_NAME,'prefix')
workDir = config.get(SECTION_NAME,'top_work_dir')
startYear = config.getint(SECTION_NAME,'start_year')
startMonth = config.getint(SECTION_NAME,'start_month')
startDay = config.getint(SECTION_NAME,'start_day')
startHour = config.getint(SECTION_NAME,'start_hour')
endYear = config.getint(SECTION_NAME,'end_year')
endMonth = config.getint(SECTION_NAME,'end_month')
endDay = config.getint(SECTION_NAME,'end_day')
endHour = config.getint(SECTION_NAME,'end_hour')
hourInterval = config.getint(SECTION_NAME,'link_interval_in_hours')

# Override from command line.
if options.prefix:
    prefix = options.prefix
    print "Overriding prefix from config file, now set as %s" %(options.prefix)

# Sanity check config settings
if not os.path.exists(workDir):
    sys.exit("FATAL, %s does not exist!" %(workDir))
if not os.access(workDir,os.R_OK | os.W_OK | os.X_OK):
    sys.exit("FATAL, %s does not have read/write/exec permission!" \
                 %(workDir))
startDateTime = datetime.datetime(startYear,startMonth,startDay,startHour)
endDateTime = datetime.datetime(endYear,endMonth,endDay,endHour)
if startDateTime > endDateTime:
    sys.exit("FATAL, start date/time is after end date/time!")
if hourInterval < 1 or hourInterval > 23:
    sys.exit("FATAL, link_interval_in_hours must be between 0 and 24!")
deltaTime = datetime.timedelta(hours=hourInterval)

# Check if starting file exists.
os.chdir(workDir)
startDateTimeString = \
    "%4.4d-%2.2d-%2.2d_%2.2d" %(startDateTime.year,
                                startDateTime.month,
                                startDateTime.day,
                                startDateTime.hour)
if not os.path.exists("%s:%s" %(prefix,startDateTimeString)):
    sys.exit("FATAL, cannot find %s WPS file for start date/time!" \
                 %(prefix))
actualDateTime = startDateTime
actualDateTimeString = startDateTimeString

# Start looping through other times
curDateTime = startDateTime + deltaTime
while curDateTime <= endDateTime:

    curDateTimeString = \
        "%4.4d-%2.2d-%2.2d_%2.2d" %(curDateTime.year,
                                    curDateTime.month,
                                    curDateTime.day,
                                    curDateTime.hour)

    # See if we have actual file for this date/time
    curFile = "%s:%s" %(prefix,curDateTimeString)
    if os.path.isfile(curFile):
        actualDateTime = curDateTime
        actualDateTimeString = curDateTimeString
        curDateTime += deltaTime
        continue

    # Delete a link if it already exists.
    if os.path.islink(curFile):
        os.remove(curFile)
    origfile = "%s:%s" %(prefix,actualDateTimeString)
    os.symlink(origfile,curFile)

    # Progress to next time        
    curDateTime += deltaTime
