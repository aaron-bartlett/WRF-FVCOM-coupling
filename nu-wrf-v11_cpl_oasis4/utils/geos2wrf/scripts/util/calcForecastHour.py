#!/usr/bin/env python

# Standard Python modules
import datetime
import sys

#------------------------------------------------------------------------------

def usage():
    txt= """
Usage: %s <startDateTime> <endDateTime>

    where each dateTime format is YYYYMMDDHH
    YYYY = year
    MM = month
    DD = day 
    HH = hour""" %(sys.argv[0])
    return

#------------------------------------------------------------------------------
def checkArguments():
    if len(sys.argv) != 3:
        usage()
        raise SystemExit, "Invalid command line arguments!"
    return sys.argv[1], sys.argv[2]

#------------------------------------------------------------------------------
def generateDateTime(timeString):
    year = int(timeString[0:4])
    month = int(timeString[4:6])
    day = int(timeString[6:8])
    hour = int(timeString[8:10])
    dt = datetime.datetime(year=year,month=month,day=day,hour=hour)
    return dt

#------------------------------------------------------------------------------
def calcForecastHours(dateTime1,dateTime2):
    delta = dateTime2 - dateTime1    
    hours= delta.days*24 + \
           delta.seconds/3600. + \
           delta.microseconds/3600000000.
    return int(hours)

#------------------------------------------------------------------------------
def main():
    timeString1,timeString2 = checkArguments()
    dt1 = generateDateTime(timeString1)
    dt2 = generateDateTime(timeString2)
    hours = calcForecastHours(dt1,dt2)
    print hours
    raise SystemExit

#------------------------------------------------------------------------------
 # Execute main only if this script is called directly.
if __name__ == "__main__":
    main()

