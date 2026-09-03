#!/usr/bin/env python
#------------------------------------------------------------------------------
# NASA/GSFC, Software Systems Support Office, Code 610.3
#------------------------------------------------------------------------------
#
# PROGRAM:  updateDateTime.py
#
# AUTHOR:
# Eric Kemp, NASA SSSO/Northrop Grumman
#
# DESCRIPTION:
# Script for incrementing a date/time by 1 day.
#
#------------------------------------------------------------------------------
"""Reads in a string from the command-line representing a date and time, and
prints a date and time incremented by 1 day to standard output.

Args:
sys.argv[1]:  String containing date/hour.  Coded as YYYYMMDDHH, where YYYY is
              the year, MM is the month, DD is the day of the month, and HH
              is the hour.

Returns:
None

Raises:
SystemExit
"""

# Standard Python modules
import datetime
import sys

#-----------------------------------------------------------------------------
def usage():
    """Prints usage of this script.

    Arg:
    None

    Returns:
    None

    Raises:
    None
    """

    txt = """
Usage: %s <currentDateTime> <incrementHours>

       where currentDateTime format is YYYYMMDDHH
       YYYY = year
       MM = month
       DD = day
       HH = hour""" %(sys.argv[0])
    return

#-----------------------------------------------------------------------------
def checkArguments():
    """Checks the command line arguments.  If an invalid number is found,
    call the usage function and terminate.  Otherwise, return the second
    command line argument, which should be the date/time string.

    Args:
    None

    Returns:
    sys.argv[1] :  Input date/time in YYYYMMDDHH format.
    sys.argv[2] :  Input hour increment.

    Raises:
    SystemExit
    """

    if len(sys.argv) != 3: 
        usage()
        raise SystemExit, "Invalid command line arguments!"
    return sys.argv[1], int(sys.argv[2])

#-----------------------------------------------------------------------------
def generateDatetime(timeString):
    """Splits the date/time string into integer year, month, day, and hour.
    Then convert to a standard datetime object and return it.

    Args:
    timeString:  A string representing date/time as YYYYMMDDHH.

    Returns:
    dt:  Python datetime object.

    Raises:
    None
    """

    year = int(timeString[0:4])
    month = int(timeString[4:6])
    day = int(timeString[6:8])
    hour = int(timeString[8:10])
    dt = datetime.datetime(year=year,month=month,day=day,hour=hour)
    return dt

#-----------------------------------------------------------------------------
def updateDatetime(dt,incrementHours):
    """Takes a Python datetime object, increments, and returns
    new object.

    Args:
    dt:  Python datetime object
    incrementHours:  Number of hours to increment

    Returns:
    New datetime object after incrementing.

    Raises:
    None
    """

    return dt + datetime.timedelta(hours=incrementHours)

#-----------------------------------------------------------------------------
def convertDatetimeToString(dt):
    """Takes a Python datetime object and constructs the stored date/time
    as a string encoded YYYYMMDDHH

    Args:
    dt:  Python datetime object

    Returns:
    timeString:  String encoded YYYYMMDDHH

    Raises:
    None
    """
    timeString = "%4.4d%2.2d%2.2d%2.2d" %(dt.year,
                                          dt.month,
                                          dt.day,
                                          dt.hour)
    return timeString

#-----------------------------------------------------------------------------
def main():
    """Main driver function.  Reads the input date/time from the command line,
    increments by 1 day, and prints updated date/time string to standard 
    output.

    Args:
    None

    Returns:
    None

    Raises:
    None
    """

    timeString,incrementHours = checkArguments()
    dt = generateDatetime(timeString)
    dt = updateDatetime(dt,incrementHours)
    timeString = convertDatetimeToString(dt)
    print timeString

#-----------------------------------------------------------------------------

# Execute main only if this script is called directly.
if __name__ == "__main__":
    main()
