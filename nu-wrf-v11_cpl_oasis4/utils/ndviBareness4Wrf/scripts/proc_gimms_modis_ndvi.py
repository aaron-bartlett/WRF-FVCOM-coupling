#!/usr/bin/env python
#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Sciences and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# SCRIPT: proc_gimms_modis_ndvi.py
#         Developed with Python 2.6.9.
#
# AUTHOR:
# Eric Kemp, NASA GSFC/SSAI
# Borrows heavily from software produced by LIS team, Code 617.
#
# DESCRIPTION:
# Script for automating processing of GIMMS MODIS-based NDVI GeoTIFF files, 
# deriving surface bareness using NDVI thresholds, and writing bareness and 
# NDVI in WPS intermediate binary format for METGRID.  Automatically downloads
# NDVI files from GIMMS FTP server unless user chooses to use local 
# pre-existing archive (e.g., the FAME directory on Discover).  Actual 
# processing of GeoTIFF data is performed by a Fortran program called by the
# script.
#
# REVISION HISTORY:
# 27 Jul 2016 - Initial version
# 03 Aug 2016 - Now supports both Aqua and Terra files.  
# 04 Aug 2016 - Improved error/exception checking.
# 01 Sep 2016 - Changed logic to use pre-staged files by default.  This makes
#               the behavior the same as the just-changed 
#               proc_sport_modis_ndvi.py script.
#
#------------------------------------------------------------------------------

# Standard Python modules
try:
    import ConfigParser
except ImportError:
    from configparser import RawConfigParser
import datetime
import ftplib
import glob
import logging
import math
import optparse
import os
import shutil
import sys

# Useful constants
SERVER_HOST = "gimms.gsfc.nasa.gov"
top_dir = {
    "TERRA" : "MODIS/std/GMOD09Q1/tif/NDVI",
    "AQUA" :  "MODIS/std/GMYD09Q1/tif/NDVI",
}
tiff_prefix = {
    "TERRA" : "GMOD09Q1",
    "AQUA"  : "GMYD09Q1",
}

SECTION_NAME = "GIMMS_MODIS"

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
# a text configuration file "gimms_modis.cfg", loops through requested days, 
# downloads and uncompresses GeoTIFF files (if necessary), and calls Fortran
# executable to process the NDVI data.
#
# This borrows heavily from the wget_gimms_modisndvi.sh script (written by
# Kristi Arsenault, Code 617) for downloading files from the GIMMS FTP server.
#
#------------------------------------------------------------------------------

def main():

    logging.info("Start NDVI processing")

    # First, check for command line flags
    logging.info('Reading command line')
    parser = optparse.OptionParser(usage="%prog [--help] [-D] [-S]")
    parser.add_option("-D","--download",action="store_true",
                      dest="downloadFiles", default=False,
                      help="Download files from remote GIMMS server")
    parser.add_option("-S","--save",action="store_true",
                      dest="saveDownloadedFiles",default=False,
                      help="Save downloaded files after processing")
    (options,args) = parser.parse_args()
    
    # Next, read config file
    logging.info('Reading config file')
    try:
        try:
            config = ConfigParser.RawConfigParser()
        except:
            config = RawConfigParser()
        config.read('gimms_modis.cfg')
        satellite = config.get(SECTION_NAME,'satellite')
        startYear = config.getint(SECTION_NAME,'start_year')
        startMonth = config.getint(SECTION_NAME,'start_month')
        startDay = config.getint(SECTION_NAME,'start_day')
        endYear = config.getint(SECTION_NAME,'end_year')
        endMonth = config.getint(SECTION_NAME,'end_month')
        endDay = config.getint(SECTION_NAME,'end_day')
        ndviBareness4WrfPath = \
            config.get(SECTION_NAME,'ndvi_bareness_4_wrf_path')
        topWorkDir = config.get(SECTION_NAME,'top_work_dir')
        if not options.downloadFiles:
            topLocalArchiveDir = \
                config.get(SECTION_NAME,'top_local_archive_dir')
        lowerLeftLat = config.getfloat(SECTION_NAME,'lower_left_lat')
        lowerLeftLon = config.getfloat(SECTION_NAME,'lower_left_lon')
        upperRightLat = config.getfloat(SECTION_NAME,'upper_right_lat')
        upperRightLon = config.getfloat(SECTION_NAME,'upper_right_lon')
        ndviBarenessThreshold1 = \
            config.getfloat(SECTION_NAME,'ndvi_bareness_threshold_1')
        ndviBarenessThreshold2 = \
            config.getfloat(SECTION_NAME,'ndvi_bareness_threshold_2')
    except ConfigParser.Error:
        logging.critical("Problem with the config file!")
        raise

    # Next, sanity check config file options
    logging.info("Sanity checking config file settings")
    satellite = satellite.upper()
    if satellite not in ["TERRA","AQUA"]:
        fatal("Invalid satellite %s selected for MODIS NDVI!" \
                  %(satellite))
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
    if startdate.year < 2000:
        fatal("Start date year is before 2000, no NDVI data available!")
    if not os.path.exists(ndviBareness4WrfPath):
        fatal("%s does not exist!" %(ndviBareness4WrfPath))
    if not os.access(ndviBareness4WrfPath,os.X_OK):
        fatal("%s is not executable!" %(ndviBareness4WrfPath))
    if not os.path.exists(topWorkDir):
        fatal("%s does not exist!" %(topWorkDir))
    if not os.access(topWorkDir,os.R_OK | os.W_OK | os.X_OK):
        fatal("%s does not have read/write/exec permission!" %(topWorkDir))
    if not options.downloadFiles:
        if not os.path.exists(topLocalArchiveDir):
            fatal("%s does not exist!" %(topLocalArchiveDir))
        if not os.access(topLocalArchiveDir,os.R_OK | os.X_OK):
            fatal("%s does not have read/exec permission!" \
                      %(topLocalArchiveDir))
    if lowerLeftLat < -60. or lowerLeftLat > 80.:
        fatal("Invalid lower left latitude!")
    if lowerLeftLon < -180. or lowerLeftLon > 180.:
        fatal("Invalid lower left longtitude!")
    if upperRightLat < -60. or upperRightLat > 80.:
        fatal("Invalid upper right latitude!")
    if upperRightLon < -180. or upperRightLon > 180.:
        fatal("Invalid upper right longtitude!")
    if lowerLeftLat > upperRightLat:
        fatal("Invalid latitude bounds!")
    if lowerLeftLon > upperRightLon:
        fatal("Invalid longitude bounds!")
    if ndviBarenessThreshold1 < 0 or ndviBarenessThreshold1 > 1:
        fatal("Invalid value for ndvi_bareness_threshold_1!")
    if ndviBarenessThreshold2 < 0 or ndviBarenessThreshold2 > 1:
        fatal("Invalid value for ndvi_bareness_threshold_2!")

    # Loop through the dates.  We'll download (if necessary) and write to WPS 
    # format.
    foundData = False
    curdate = startdate    
    while curdate <= enddate:
        year = curdate.year
        jday = curdate.timetuple()[7]

        # Search for pre-existing local files, unless downloading is
        # specified at command line.
        if not options.downloadFiles:
            tmpArchiveDir = "%s/%4.4d" %(topLocalArchiveDir,year)
            if not os.path.exists(tmpArchiveDir):
                fatal("Cannot access %s" %(tmpArchiveDir))
            tmpArchiveDir = "%s/%4.4d/%3.3d" %(topLocalArchiveDir,year,jday)
            if not os.path.exists(tmpArchiveDir):
                logging.warning("No %s MODIS NDVI data available for %s" \
                                    %(satellite,curdate))
                curdate += datetime.timedelta(days=1)        
                continue
            files = glob.glob("%s/%s*" %(tmpArchiveDir,tiff_prefix[satellite]))
            if len(files) == 0:
                logging.warning("No %s MODIS NDVI data available for %s" \
                                    %(satellite,curdate))
                curdate += datetime.timedelta(days=1)
                continue                              
            foundData = True

        else:

            # Open anonymous FTP connection
            try:
                ftp = ftplib.FTP(SERVER_HOST)
                ftp.login() 
            except ftplib.all_errors:
                logging.critical("Problem logging into GIMMS FTP server!")
                raise

            # See if requested year is available
            dirname = top_dir[satellite] + "/%4.4d" %(year)
            try:
                ftp.cwd(dirname)
            except ftplib.all_errors:
                logging.critical( \
                    "Cannot access directory on GIMMS FTP server!" %(dirname))
                raise

            # See if requested day of year is available
            try:
                jdays = ftp.nlst()
            except ftplib.all_errors:
                logging.critical( \
                    "Cannot get directory contents on GIMMS FTP server!")
                raise
            s_jday = "%3.3d" %(jday)
            if s_jday not in jdays:
                logging.warning("No NDVI data available for %s" %(curdate))
                try:
                    ftp.quit()
                except ftplib.all_errors:
                    logging.critical( \
                        "Cannot close connection to GIMMS FTP server!")
                    raise
                curdate += datetime.timedelta(days=1)        
                continue
            dirname = "%3.3d" %(jday)
            try:
                ftp.cwd(dirname)
            except ftplib.all_errors:
                logging.critical(\
                    "Cannot access directory on GIMMS FTP server!")
                raise

            # Create local output directory
            tmpDownloadDir = "%s/download/%4.4d/%3.3d" %(topWorkDir,year,jday)
            if not os.path.exists(tmpDownloadDir):
                try:
                    os.makedirs(tmpDownloadDir)
                except os.error:
                    logging.critical("Problem creating directory %s" \
                                         %(tmpDownloadDir))
                    raise
            try:
                os.chdir(tmpDownloadDir)
            except os.error:
                logging.critical("Problem accessing directory %s" \
                                     %(tmpDownloadDir))
                raise

            # Loop through the files and download the GZIPPED GeoTIFF files.
            # Only select the tiles needed to fill the requested lat/lon 
            # bounds.
            try:
                remotefiles_all = ftp.nlst()
            except ftplib.all_errors:
                logging.critical("Problem listing files on GIMMS FTP server!")
                raise
            remotefiles = []
            prefix = "%s.A%4.4d%3.3d.08d.latlon." \
                %(tiff_prefix[satellite],year,jday)
            west_x  = int(math.floor((180 + lowerLeftLon ) / 9))
            south_y = int(math.floor(( 90 - lowerLeftLat ) / 9))
            east_x  = int(math.floor((180 + upperRightLon) / 9))
            north_y = int(math.floor(( 90 - upperRightLat) / 9))
            for x in range(west_x,east_x+1):
                for y in range(north_y,south_y+1):
                    remotefile = prefix + \
                        "x%2.2dy%2.2d.6v1.NDVI.tif.gz" %(x,y)
                    if remotefile not in remotefiles_all:
                        logging.warning("%s is not on remote server" \
                                            %(remotefile))
                        continue
                    remotefiles.append(remotefile)
                    try:
                        fhandle = open(remotefile,'wb')
                    except IOError:
                        logging.critical("Cannot open %s!" %(remotefile))
                        raise
                    logging.info("Fetching %s" %(remotefile))
                    try:
                        ftp.retrbinary('RETR %s' %(remotefile),fhandle.write)
                    except ftplib.all_errors:
                        logging.critical("Cannot download %s!" %(remotefile))
                        raise
                    fhandle.close()
                     
            # Close the FTP connection.
            try:
                ftp.quit()
            except ftplib.all_errors:
                logging.critical(\
                    "Problem closing connection to GIMMS FTP server")
                raise

            # Clean up for this day.
            del remotefiles_all
            if len(remotefiles) == 0:
                logging.warning("No NDVI files found for %s " %(curdate))
                curdate += datetime.timedelta(days=1)        
                continue
            del remotefiles
            foundData = True

            # Now uncompress the files
            logging.info("Running gunzip...")
            cmd = 'gunzip *.tif.gz'
            try:
                status = os.system(cmd)
                if status != 0:
                    fatal("Problem running gunzip on GeoTIFF files!")
            except os.error:
                logging.critical("Problem running gunzip on GeoTIFF files!")
                raise

        # Now create namelist file for ndviBareness4Wrf
        try:
            os.chdir(topWorkDir)
        except os.error:
            logging.critical("Cannot access work directory %s!" %(topWorkDir))
            raise
        logging.info("Creating namelist.ndviBareness4Wrf")
        try:
            fd = open("namelist.ndviBareness4Wrf","w")
            fd.write("\n &settings\n")
            if satellite == "AQUA":
                fd.write("   ndviFormat = 10\n")
            elif satellite == "TERRA":
                fd.write("   ndviFormat = 11\n")
            fd.write("   lowerLeftLat = %f\n" %(lowerLeftLat))
            fd.write("   lowerLeftLon = %f\n" %(lowerLeftLon))
            fd.write("   upperRightLat = %f\n" %(upperRightLat))
            fd.write("   upperRightLon = %f\n" %(upperRightLon))
            fd.write("   ndviBarenessThreshold1 = %f\n" \
                         %(ndviBarenessThreshold1))
            fd.write("   ndviBarenessThreshold2 = %f\n" \
                         %(ndviBarenessThreshold2))
            if not options.downloadFiles:
                fd.write("   inputDirectory = '%s'\n" %(topLocalArchiveDir))
            else:
                fd.write("   inputDirectory = '%s/download'\n" %(topWorkDir))
            fd.write("   year = %d\n" %(curdate.year))
            fd.write("   month = %d\n" %(curdate.month))
            fd.write("   dayOfMonth = %d\n" %(curdate.day))
            day_of_year = datetime.date(curdate.year,curdate.month,curdate.day)
            doy = day_of_year.strftime('%j')
            fd.write("   dayOfYear = %s\n" %(doy))
            fd.write("   outputDirectory = './'\n")
            fd.write("   outputFilePrefix = 'GIMMS_MODIS'\n")
            fd.write("   /\n")
            fd.close()
        except IOError:
            logging.critical("Problem writing to namelist.ndviBareness4Wrf!")
            raise

        # Now run ndviBareness4Wrf
        logging.info("Calling ndviBareness4Wrf...")
        try:
            status = os.system(ndviBareness4WrfPath)
            if status != 0:
                fatal("Problem running ndviBareness4Wrf!")
        except os.error:
            logging.critical("Problem running ndviBareness4Wrf!")
            raise

        # Optionally clean up downloaded files.  Pre-existing files will not
        # be affected.
        if options.downloadFiles:
            if not options.saveDownloadedFiles:
                logging.debug("Removing %s/download" %(topWorkDir))
                try:
                    shutil.rmtree("%s/download" %(topWorkDir))
                except shutil.Error:
                    logging.critical("Problem deleting %s/download!" \
                                         %(topWorkDir))
                    raise
            else:
                logging.debug("Saving %s" %(tmpDownloadDir))

        # Now move to next day
        curdate += datetime.timedelta(days=1)        

    # The end
    logging.info("End NDVI processing")
    if not foundData:
        fatal("No %s MODIS NDVI files were processed!" %(satellite))

# Only invoke main() if this script is called directly.
if __name__ == "__main__":
    main()
    
