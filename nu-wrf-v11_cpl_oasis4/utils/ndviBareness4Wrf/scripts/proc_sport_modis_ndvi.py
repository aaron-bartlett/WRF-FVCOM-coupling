#!/usr/bin/env python
#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Sciences and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# SCRIPT: proc_sport_modis_ndvi.py
#         Developed with Python 2.6.9.
#
# AUTHOR:
# Eric Kemp, NASA GSFC/SSAI
#
# DESCRIPTION:  Automates processing of SPORT MODIS-based NDVI GrADS files,
# deriving surface bareness using NDVI thresholds, and writing bareness and
# NDVI in WPS intermediate binary format for METGRID.  Actual processing of
# GrADS data is performed by a Fortran program called by the script.
#
# REVISION HISTORY:
# 29 Jul 2016 - Initial version.
# 04 Aug 2016 - Now reports error if no files are found.  Improved error/
#               exception checking.
# 01 Sep 2016 - Changed code to always used pre-fetched files (no download
#               server exists and probably never will).  Fixed typo with
#               checking for exceptions from ConfigParser.
#
#------------------------------------------------------------------------------

# Standard python modules
import ConfigParser
import datetime
import logging
import os
import shutil
import sys

# Useful constants
SECTION_NAME = "SPORT_MODIS"

# Configure logging for the script.
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', \
                        level=logging.INFO)

# Function for handling fatal errors
def fatal(message):
    logging.critical(message)
    sys.exit(1)

#------------------------------------------------------------------------------
#
# FUNCTION: main
#
# DESCRIPTION:  Main driver for this script.  Processes command line arguments,
# a text configuration file "sport_modis.cfg", loops through requested days, 
# and calls Fortran executable to process the NDVI data.
#
#------------------------------------------------------------------------------

def main():
    
    logging.info("Start NDVI processing")

    # Read config file
    logging.debug('Reading config file')
    try:
        config = ConfigParser.RawConfigParser()
        config.read('sport_modis.cfg')
        startYear = config.getint(SECTION_NAME,'start_year')
        startMonth = config.getint(SECTION_NAME,'start_month')
        startDay = config.getint(SECTION_NAME,'start_day')
        endYear = config.getint(SECTION_NAME,'end_year')
        endMonth = config.getint(SECTION_NAME,'end_month')
        endDay = config.getint(SECTION_NAME,'end_day')
        ndviBareness4WrfPath = config.get(SECTION_NAME, \
                                              'ndvi_bareness_4_wrf_path')
        topWorkDir = config.get(SECTION_NAME,'top_work_dir')
        topLocalArchiveDir = config.get(SECTION_NAME, \
                                            'top_local_archive_dir')
        ndviBarenessThreshold1 = \
            config.getfloat(SECTION_NAME,'ndvi_bareness_threshold_1')
        ndviBarenessThreshold2 = \
            config.getfloat(SECTION_NAME,'ndvi_bareness_threshold_2')
        outputFilePrefix = config.get(SECTION_NAME,'output_file_prefix')
    except ConfigParser.Error:
        logging.critical("Problem with the config file!")
        raise

    # Next, sanity check config file options
    logging.debug('Sanity checking config file settings')
    try:
        startdate = datetime.date(startYear,startMonth,startDay)
    except ValueError:
        logging.critical("Invalid start date!")
        raise
    try:
        enddate = datetime.date(endYear,endMonth,endDay)
    except ValueError:
        logging.critical("Invalid end date!")
        raise
    if enddate < startdate:
        fatal("Start date is after end date!")
    if startdate.year < 2010:
        fatal("Start date year is before 2010, no NDVI data available!")
    if not os.path.exists(ndviBareness4WrfPath):
        fatal("%s does not exist!" %(ndviBareness4WrfPath))
    if not os.access(ndviBareness4WrfPath,os.X_OK):
        fatal("%s is not executable!" %(ndviBareness4WrfPath))
    if not os.path.exists(topWorkDir):
        fatal("%s does not exist!" %(topWorkDir))
    if not os.access(topWorkDir,os.R_OK | os.W_OK | os.X_OK):
        fatal("%s does not have read/write/exec permission!" %(topWorkDir))
    if not os.path.exists(topLocalArchiveDir):
        fatal("%s does not exist!" %(topLocalArchiveDir))
    if not os.access(topLocalArchiveDir,os.R_OK | os.X_OK):
        fatal("%s does not have read/exec permission!" \
                  %(topLocalArchiveDir))
    
    # Loop through the dates.  
    foundData = False
    curdate = startdate
    while curdate <= enddate:
        year = curdate.year
        month = curdate.month
        dayOfMonth = curdate.day

        tmpArchiveDir = "%s/%4.4d" %(topLocalArchiveDir,year)
        if not os.path.exists(tmpArchiveDir):
            fatal("Cannot access %s" %(tmpArchiveDir))
        tmpArchiveDir = "%s/%4.4d/%2.2d" %(topLocalArchiveDir,year,month)
        if not os.path.exists(tmpArchiveDir):
            fatal("Cannot access %s" %(tmpArchiveDir))
            
        # Make sure local file exists.  If necessary, use gunzip.
        localFile = "%s/ndvi1km_lismask_%4.4d%2.2d%2.2d.bin" \
            %(tmpArchiveDir,year,month,dayOfMonth)
        if not os.path.exists(localFile):
            localFileGz = localFile + ".gz"
            if not os.path.exists(localFileGz):
                logging.warning("No NDVI data available for %s" %(curdate))
                continue
            try:
                os.chdir(tmpArchiveDir)
            except os.error:
                logging.critical("Cannot access %s!" %(tmpArchiveDir))
                raise
            cmd = 'gunzip %s' %(localFileGz)
            logging.debug('Running %s' %(cmd))
            try:
                status = os.system(cmd)
                if status != 0:
                    fatal("Problem running gunzip!")
            except os.error:
                logging.critical("Problem running gunzip!")
                raise
        foundData = True
            
        # Now create namelist file for ndviBareness4Wrf
        try:
            os.chdir(topWorkDir)
        except os.error:
            logging.critical("Cannot access work directory %s!" %(topWorkDir))
            raise
        logging.debug("Creating namelist.ndviBareness4Wrf")
        try:
            fd = open("namelist.ndviBareness4Wrf","w")
            fd.write("\n &settings\n")
            fd.write("   ndviFormat = 1\n")
            fd.write("   ndviBarenessThreshold1 = %f\n" \
                         %(ndviBarenessThreshold1))
            fd.write("   ndviBarenessThreshold2 = %f\n" \
                         %(ndviBarenessThreshold2))
            fd.write("   inputDirectory = '%s'\n" %(topLocalArchiveDir))
            fd.write("   year = %d\n" %(curdate.year))
            fd.write("   month = %d\n" %(curdate.month))
            fd.write("   dayOfMonth = %d\n" %(curdate.day))
            fd.write("   outputDirectory = './'\n")
            fd.write("   outputFilePrefix = '%s'\n" %(outputFilePrefix))
            fd.write("   /\n")
            fd.close()
        except IOError:
            logging.critical("Problem writing to namelist.ndviBareness4Wrf")
            raise

        # Now run ndviBareness4Wrf
        logging.debug("Calling ndviBareness4Wrf...")
        try:
            status = os.system(ndviBareness4WrfPath)
            if status != 0:
                fatal("Problem running ndviBareness4Wrf!")
        except os.error:
            logging.critical("Problem running ndviBareness4Wrf!")
            raise

        # Now move to next day
        curdate += datetime.timedelta(days=1)

    # The end
    logging.info("End NDVI processing")
    if not foundData:
        fatal("No SPoRT MODIS NDVI files were processed!" \
                  %(satellite))

# Only invoke main() if this script is called directly.
if __name__ == "__main__":
    main()

