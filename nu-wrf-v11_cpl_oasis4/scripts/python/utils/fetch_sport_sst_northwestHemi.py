#!/usr/bin/env python
# ------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Sciences and Technology Office,
# Code 606
# ------------------------------------------------------------------------------
#
# SCRIPT: fetch_sport_sst_northwestHemi.py
#         Developed with Python 2.6.9.
#
# AUTHOR:  Eric Kemp, NASA GSFC/SSAI
#
# DESCRIPTION:  Script for downloading NASA SPoRT SST GRIB2 data for NorthWest
# hemispheric region.
#
# REVISION HISTORY:  31 Aug 2016 - Initial version
#
# ------------------------------------------------------------------------------

# Standard Python modules
import configparser
import datetime
import ftplib
import logging
import optparse
import os
import sys

# Useful constants
SERVER_HOST = "geo.msfc.nasa.gov"
SST_DIR = "SPoRT/sst/northwestHemi/grib2"
SECTION_NAME = "SPORT_SST"

# Configure logging for the script.  Log entries should include the local
# time on the machine.  Messages with severity below INFO will not be written.
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s", level=logging.INFO
)

# Function for handling fatal errors
def fatal(message):
    logging.critical(message)
    sys.exit(1)


# ------------------------------------------------------------------------------
#
# FUNCTION: main
#
# DESCRIPTION:  Main driver for this script.  Processes command line arguments,
# a text configuration file "sport_sst_northwestHeni.cfg", loops through
# requested days and times, and downloads and uncompresses GRIB2 files.
#
# ------------------------------------------------------------------------------


def main():

    logging.info("Start SPORT SST download")

    # First, check for command line flags
    logging.debug("Reading command line")
    parser = optparse.OptionParser(usage="%prog [--help] [-L]")
    parser.add_option(
        "-L",
        "--latest",
        action="store_true",
        dest="fetchLatestFiles",
        default=False,
        help="Override cfg file, fetch prior 24 hours of files",
    )
    (options, args) = parser.parse_args()

    # Next, read config file
    logging.debug("Reading config file")
    try:
        config = configparser.ConfigParser.RawConfigParser()
        config.read("sport_sst_northwestHemi.cfg")
        topDownloadDir = config.get(SECTION_NAME, "top_download_dir")
    except configparser.ConfigParser.Error:
        logging.critical("Problem with the config file!")
        raise
    if not options.fetchLatestFiles:
        try:
            startYear = config.getint(SECTION_NAME, "start_year")
            startMonth = config.getint(SECTION_NAME, "start_month")
            startDay = config.getint(SECTION_NAME, "start_day")
            startHour = config.getint(SECTION_NAME, "start_hour")
            endYear = config.getint(SECTION_NAME, "end_year")
            endMonth = config.getint(SECTION_NAME, "end_month")
            endDay = config.getint(SECTION_NAME, "end_day")
            endHour = config.getint(SECTION_NAME, "end_hour")
        except configparser.ConfigParser.Error:
            logging.critical("Problem with the config file!")
            raise

    # Next, sanity check config file options
    logging.debug("Sanity checking config file settings")
    if not os.path.exists(topDownloadDir):
        fatal("%s does not exist!" % (topDownloadDir))
    if not os.access(topDownloadDir, os.R_OK | os.W_OK | os.X_OK):
        fatal("%s does not have read/write/exec permission!" % (topDownloadDir))
    if not options.fetchLatestFiles:
        try:
            startdatetime = datetime.datetime(
                startYear, startMonth, startDay, startHour
            )
        except ValueError:
            logging.critical("Invalid start datetime!")
            raise
        try:
            enddatetime = datetime.datetime(endYear, endMonth, endDay, endHour)
        except ValueError:
            logging.critical("Invalid end datetime!")
            raise
        if enddatetime < startdatetime:
            fatal("Start datetime is after end datetime!")
        if startdatetime.hour not in [6, 18]:
            fatal("Start hour must be 6 or 18 UTC!")
        if enddatetime.hour not in [6, 18]:
            fatal("End hour must be 6 or 18 UTC!")

    # If fetching latest files, dynamically set start and end date/times.
    # The end date/time will be the more recent previous SST time, and
    # the start date/time will be 24 hours prior.
    if options.fetchLatestFiles:
        logging.info("Resetting download dates for latest data")
        enddatetime = datetime.datetime.utcnow()
        enddatetime = enddatetime.replace(minute=0, second=0, microsecond=0)
        if enddatetime.hour > 18:
            enddatetime = enddatetime.replace(hour=18)
        elif enddatetime.hour > 6:
            enddatetime = enddatetime.replace(hour=6)
        elif enddatetime.hour < 6:
            # Set to 18 UTC for the previous day
            enddatetime = enddatetime - datetime.timedelta(days=1)
            enddatetime = enddatetime.replace(hour=18)
        startdatetime = enddatetime - datetime.timedelta(days=1)

    # Open anonymous FTP connection
    try:
        ftp = ftplib.FTP(SERVER_HOST)
        ftp.login()
    except ftplib.all_errors:
        logging.critical("Problem logging into SPoRT FTP server!")
        raise

    # Go to appropriate directory
    try:
        ftp.cwd(SST_DIR)
    except ftplib.all_errors:
        logging.critical("Cannot access directory on SPoRT FTP server!")
        raise

    # Get list of GRIB2 files
    try:
        remotefiles_all = ftp.nlst()
    except ftplib.all_errors:
        logging.critical("Problem listing files on SPoRT FTP server!")
        raise
    remotefiles = [file for file in remotefiles_all if "sstcomp" in file]
    if len(remotefiles) == 0:
        fatal("No SST files found on SPoRT FTP server!")

    # Loop through the dates and times to download.
    foundData = False
    curdatetime = startdatetime
    deltaHours = datetime.timedelta(hours=12)
    while curdatetime <= enddatetime:

        # Is file for current date/time on server?
        year = curdatetime.year
        month = curdatetime.month
        day = curdatetime.day
        hour = curdatetime.hour
        curfile = "%4.4d%2.2d%2.2d_%2.2d00_sport_nhemis_sstcomp.grb2.gz" % (
            year,
            month,
            day,
            hour,
        )
        if curfile not in remotefiles:
            logging.warning("Cannot find %s on SPoRT FTP server!" % (curfile))
            curdatetime += deltaHours
            continue

        # Move to download directory.
        tmpDownloadDir = "%s/Y%4.4d/M%2.2d" % (
            topDownloadDir,
            curdatetime.year,
            curdatetime.month,
        )
        if not os.path.exists(tmpDownloadDir):
            try:
                os.makedirs(tmpDownloadDir)
            except os.error:
                logging.critical("Problem accessing directory %s" % (tmpDownloadDir))
                raise
        os.chdir(tmpDownloadDir)

        # Skip download if file is already present.
        localFile = curfile.split(".gz")[0]
        if os.path.exists(localFile):
            logging.info(
                "GRIB2 file %s already in local download directory" % (localFile)
            )
            curdatetime += deltaHours
            continue

        # Download file
        try:
            fhandle = open(curfile, "wb")
        except IOError:
            logging.critical("Cannot open %s!" % (curfile))
            raise
        logging.info("Fetching %s" % (curfile))
        try:
            ftp.retrbinary("RETR %s" % (curfile), fhandle.write)
        except ftplib.all_errors:
            logging.critical("Cannot download %s!" % (curfile))
            raise
        fhandle.close()

        # Uncompress the new file
        logging.info("Running gunzip...")
        cmd = "gunzip %s" % (curfile)
        try:
            status = os.system(cmd)
            if status != 0:
                fatal("Problem running gunzip on %s" % (curfile))
        except os.error:
            logging.critical("Problem running gunzip on %s" % (curfile))
            raise

        # Next date/time
        foundData = True
        curdatetime += deltaHours

    # Close the FTP connection
    try:
        ftp.quit()
    except ftplib.all_errors:
        logging.critical("Problem closing connection to SPoRT FTP server")
        raise

    # The end
    logging.info("End SPoRT SST download")
    if not foundData:
        fatal("No SPoRT SST files were downloaded!")


# Only invoke main() if this script is called directly.
if __name__ == "__main__":
    main()
