#!/usr/bin/env python
#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Sciences and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# SCRIPT: link_gimms_wps_files.py
#
# AUTHOR:
# Eric Kemp, NASA GSFC/SSAI
#
# DESCRIPTION:
# Creates symbolic links to WPS intermediate files containing NASA GIMMS 
# MODIS-based NDVI and bareness.  Although the WPS format includes internal
# time records for individual fields, METGRID does not use them, and instead
# assumes the valid date/time matches the date/time used in the file name.
# Thus symbolic links can be used to force METGRID to include these fields
# in time-varying lateral boundary conditions without wasting disk space.
#
# REVISION HISTORY:
# 29 Jul 2016 - Initial version
# 03 Aug 2016 - Rewrote code and renamed script.  Time interval for creating 
#               links is now inclusive.  More careful checking is made to 
#               discriminate between existing files and existing links.  If 
#               valid files are found while looping through the full time 
#               period, these files will become new targets for linking. Also
#               discriminates between Aqua and Terra MODIS files.
#
#------------------------------------------------------------------------------

# Standard Python modules
try:
    import ConfigParser
except ImportError:
    from configparser import RawConfigParser
import datetime
import glob
import os
import sys

# Useful constants
SECTION_NAME = "LINK_GIMMS"

# Read config file
try:
    config = ConfigParser.RawConfigParser()
except:
    config = RawConfigParser()
config.read("link_gimms.cfg")
satellite = config.get(SECTION_NAME,'satellite')
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

# Sanity check config settings
satellite = satellite.upper()
if satellite not in ["TERRA","AQUA"]:
    sys.exit("FATAL, invalid satellite %s selected for MODIS NDVI!" \
                 %(satellite))
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

# Get list of prefixes for existing files (one prefix per tile)
os.chdir(workDir)
startDateTimeString = \
    "%4.4d-%2.2d-%2.2d_%2.2d" %(startDateTime.year,
                                startDateTime.month,
                                startDateTime.day,
                                startDateTime.hour)
startFiles = glob.glob("%s*:%s" %(satellite,startDateTimeString))
prefixes = []
for startFile in startFiles:
    if not os.path.isfile(startFile):
        continue
    if os.path.islink(startFile):
        continue
    prefixes.append(startFile.split(":")[0])
if len(prefixes) == 0:
    sys.exit("FATAL, no GIMMS %s WPS files found for start date/time!" \
                 %(satellite))
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

    # See if we have actual files for this date/time
    curFiles = glob.glob("%s*:%s" %(satellite,curDateTimeString))
    curPrefixes = []
    for curFile in curFiles:
        if not os.path.isfile(curFile):
            continue
        if os.path.islink(curFile):
            continue
        curPrefixes.append(curFile.split(":")[0])
    if len(curPrefixes) != 0:
        prefixes = curPrefixes
        actualDateTime = curDateTime
        actualDateTimeString = curDateTimeString
        curDateTime += deltaTime
        continue

    # No actual files found in current date/time, so create the links.
    for prefix in prefixes:
        newlink = "%s:%s" %(prefix,curDateTimeString)
        # Delete a link if it already exists.
        if os.path.islink(newlink):
            os.remove(newlink)
        origfile = "%s:%s" %(prefix,actualDateTimeString)
        os.symlink(origfile,newlink)        

    # Progress to next time        
    curDateTime += deltaTime

