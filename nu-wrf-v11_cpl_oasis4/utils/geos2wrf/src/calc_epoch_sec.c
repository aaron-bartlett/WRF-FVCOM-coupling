//-----------------------------------------------------------------------------
// NASA/GSFC, Computational and Information Science and Technology Office,
// Code 606
//-----------------------------------------------------------------------------
//
// MODULE: calc_epoch_sec
//
// AUTHOR: Eric Kemp, NASA CISTO/SSAI
//
// DESCRIPTION: Calculates epoch time from given date and time in UTC. Acts
// as wrapper function for calling C standard libary from Fortran.
//
// REVISION:
// 20 May 2015 - First version
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Standard C library headers
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//-----------------------------------------------------------------------------
// FUNCTION: calc_epoch_sec
//
// DESCRIPTION: Calculates epoch time from given date and time in UTC. Acts 
// as wrapper function for calling C standard libary from Fortran.
//-----------------------------------------------------------------------------

int calc_epoch_sec(int * year, int * month, int * day,
		   int * hour, int * minute, int * second) {

  // Local variables
  struct tm timeObj;
  time_t epochSec;

  // Populate the date and time information of a struct tm type, then
  // pass to mktime to get the epoch seconds.
  timeObj.tm_year = (*year) - 1900;
  timeObj.tm_mon  = (*month) - 1;
  timeObj.tm_mday = *day;
  timeObj.tm_hour = *hour;
  timeObj.tm_min  = *minute;
  timeObj.tm_sec  = *second;

  timeObj.tm_isdst = 0; //No daylight savings time.

  epochSec = mktime(&timeObj);
  if (epochSec == -1) { 
    printf("ERROR calculating epoch seconds from input date/time!\n");
    exit(EXIT_FAILURE);
  }

  // time_t is C implementation dependent. Convert to more standard int before
  // passing back to Fortran.
  return (int) epochSec;

}
