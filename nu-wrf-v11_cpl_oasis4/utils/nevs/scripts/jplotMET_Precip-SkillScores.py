#!/usr/bin/env python
#------------------------------------------------------------------------------
# AUTHORS:      Bruce Van Aartsen
# AFFILIATION:  NASA GSFC / SSAI
# DATE:         May 27, 2016
#
# DESCRIPTION:
# This script will parse *.cts output files from the MET GridStat script to 
# graph a time-series of 3 Skill Scores  for given parameters
# 
# This script is meant to be executed by Jupyter Notebook: NUWRF_validation.ipynb
# It requires the following parameters to be predefined in the Notebook environment:
# fileDate - Python list of analysis periods in this format: YYYYMMDDHH
# inputRootDir - Directory that contains the MET PointStat output *.cnt files
# nuwrfRunDirs - Python list of subdirectories containing NU-WRF ensemble runs
# nuwrfRunLabels - Python list of labels corresponding to above ensemble runs
#                                   
#------------------------------------------------------------------------------

import os
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter
from matplotlib.dates import DateFormatter, drange, DayLocator
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pylab
import datetime
import csv

#------------------------------------------------------------------------------
# Set up date/time variables for the experiment 
#------------------------------------------------------------------------------

bDate = fileDate[1]
beginYear = bDate[:4]
beginMon = bDate[4:6]
beginDate = bDate[6:8]
beginHour = bDate[-2:]
eDate = fileDate[-1]
endYear = eDate[:4]
endMon = eDate[4:6]
endDate = eDate[6:8]
endHour = eDate[-2:]
hourIncr = 24

#------------------------------------------------------------------------------
# Set up string variables for the input files/data
#------------------------------------------------------------------------------

startFilename = "grid_stat_APCP_24_"
midFilename = "0000L_"
endFilename = "0000V_cts.txt"
#inputRootDir = "/discover/nobackup/projects/nu-wrf/members/bvanaart/SkillScoreRuns/MET_output/GridStat/round"
#exptDirList = ["0_default","2_bl_pbl_1_sf_sfclay_11","2_bl_pbl_2_sf_sfclay_2","3_ra_lw_phys_4_sw_phys_4","3_ra_lw_phys_56_sw_phys_56","4_moist_adv_opt_2","4_moist_adv_opt_4"]
exptDirList = nuwrfRunDirs

#exptLabels = [" Def_Goddard4"," YSU PBL, MM5 M-O SL", " MYJ TKE PBL, M-O-J SL", " RRTMG"," 2014 Goddard rad"," Monotonic advec"," 5th-order WENO advec"]
exptLabels = nuwrfRunLabels
myColors = [colorConverter.to_rgba(c)
    for c in ('red','black','goldenrod','blue','darkseagreen','darkgreen','mediumpurple','indigo','darkturquoise')]
numExpts = len(exptLabels)

#Set output PDF file name, which will host several plots
eLabel = '-'.join(str(x) for x in exptLabels[0:2])
pdfName = 'Precip-SkillScores_' + ''.join(eLabel.split()) + '.pdf'
pp = PdfPages(pdfName)

#------------------------------------------------------------------------------
# Find total # of time-periods to plot
#------------------------------------------------------------------------------

startDate = datetime.datetime(int(beginYear), int(beginMon), int(beginDate), int(beginHour))
endDate = datetime.datetime(int(endYear), int(endMon), int(endDate), int(endHour))
deltaTimePeriod = endDate - startDate
#print deltaTimePeriod

deltaHours = deltaTimePeriod.total_seconds() / 3600
numPeriods = int((deltaHours / int(hourIncr)) + 1)
print ('numPeriods = %d' %(numPeriods))

#Create list of fcst valid times
dateIncr = datetime.timedelta(hours = int(hourIncr))
print (dateIncr)
allDates = np.array([startDate + datetime.timedelta(hours = int(hourIncr*i)) for i in range(numPeriods)])
#Double check:
for dateTime in allDates:
    print (dateTime.strftime("%Y-%m-%d-%HZ"))

#------------------------------------------------------------------------------
# Set up arrays to hold stats values
#------------------------------------------------------------------------------

GSSPrecip = np.zeros((numExpts,numPeriods))
HKPrecip = np.zeros((numExpts,numPeriods))
HSSPrecip = np.zeros((numExpts,numPeriods))
GSSPrecip50 = np.zeros((numExpts,numPeriods))
HKPrecip50 = np.zeros((numExpts,numPeriods))
HSSPrecip50 = np.zeros((numExpts,numPeriods))

#------------------------------------------------------------------------------
# Begin looping over the experiments and time periods
# *** Note: Initial fcstHH is assumed to be equal to hourIncr (not 0) ***
#------------------------------------------------------------------------------

for exptNum in range(0, numExpts):
   for periodNum in range(0, numPeriods):
      date = startDate + datetime.timedelta(hours=(hourIncr*periodNum))
      fcstHH = str((periodNum + 1) * hourIncr)
      thisYY = date.strftime('%Y')
      thisMM = date.strftime('%m')
      thisDD = date.strftime('%d')
      thisHH = date.strftime('%H')
      pathName = inputRootDir + exptDirList[exptNum] + "/" + startFilename + fcstHH + midFilename + thisYY + thisMM + thisDD + "_" + thisHH + endFilename
      #print 'path = ' + pathName

#------------------------------------------------------------------------------
# Parse *cts file to extract desired stats. Most useful:
# Column#  Name         Description
# -------  ----         -----------
# 9        FCST_VAR     Model variable (i.e., T2, Q2, U10, V10)
# 15       INTERP_MTHD  Interpolation method (NEAREST, MEDIAN, DW_MEAN)
# 17       FCST_THRESH  The threshold applied to the forecast
# 22       TOTAL        Total number of matched pairs for comparison
# 66       GSS          Gilbert Skill Score
# 69       HK           Hanssen-Kuipers Discriminant (True skill statistic)
# 74       HSS          Heidke Skill Score
#------------------------------------------------------------------------------

      file = open(pathName, 'r')
      # Use csv.reader to read each line and split it into columns
      for line in csv.reader(file, delimiter=" ", skipinitialspace=True):
         # Check column 9 for RAINNC0_24 and col 14 for VX_MASK=WRF_En (use column# - 1)s
         if (line[8] == "RAINNC0_24") and (line[13] == "WRF_Ens"):
            
            # Get Skill Scores from lines where FCST_THRESH >10 
            # Check column 17 for FCST_THRESH
            if (line[16] == ">10.0"):
               GSSPrecip[exptNum,periodNum] = float(line[65])
               HKPrecip[exptNum,periodNum] = float(line[68])
               HSSPrecip[exptNum,periodNum] = float(line[73])
               
            # Get Skill Scores from lines where FCST_THRESH >50 
            # Check column 17 for FCST_THRESH
            if (line[16] == ">50.0"):
               GSSPrecip50[exptNum,periodNum] = float(line[65])
               HKPrecip50[exptNum,periodNum] = float(line[68])
               HSSPrecip50[exptNum,periodNum] = float(line[73])
            
#print GSSPrecip

#------------------------------------------------------------------------------
# Start plotting the data
#------------------------------------------------------------------------------

#Create doc to host several plots (now defined near top)
#pp = PdfPages('Precip-SkillScores.pdf')

#GSS plots - Precip > 10mm
fig = plt.figure()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, GSSPrecip[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('Gilbert Skill Score - Precip')
plt.title('Gilbert Skill Score - 24-hr Precip > 10mm ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
# save to PDF document (defined near the top)
pp.savefig()
pylab.savefig('Precip_GilbertSkillScore.png')

#HK plots - Precip > 10mm
#fig = plt.figure()
plt.clf()

for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HKPrecip[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('HK Discriminant - 24-hr Precip')
plt.title('Hanssen-Kuipers (True Skill) - 24-hr Precip > 10mm ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Precip_H-KDiscriminant.png')

#HSS plots - Precip > 10mm
#fig = plt.figure()
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HSSPrecip[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('Heidke Skill Score - Precip')
plt.title('Heidke Skill Score - 24-hr Precip > 10mm ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('Precip_HeidkeSkillScore.png')

#GSS plots - Precip > 50mm
#fig = plt.figure()
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, GSSPrecip50[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('Gilbert Skill Score - Precip')
plt.title('Gilbert Skill Score - 24-hr Precip > 50mm ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
# save to PDF document (defined near the top)
pp.savefig()
pylab.savefig('Precip50mm_GilbertSkillScore.png')

#HK plots - Precip > 50mm
#fig = plt.figure()
plt.clf()

for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HKPrecip50[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('HK Discriminant - 24-hr Precip')
plt.title('Hanssen-Kuipers (True Skill) - 24-hr Precip > 50mm ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Precip50mm_H-KDiscriminant.png')


#HSS plots - Precip > 50mm
#fig = plt.figure()
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HSSPrecip50[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('Heidke Skill Score - Precip')
plt.title('Heidke Skill Score - 24-hr Precip > 50mm ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()


#Close PDF doc
pp.close()
pylab.savefig('Precip50mm_HeidkeSkillScore.png')
plt.close()
