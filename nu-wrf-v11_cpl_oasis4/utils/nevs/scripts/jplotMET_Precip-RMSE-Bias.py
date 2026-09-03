#!/usr/bin/env python
#------------------------------------------------------------------------------
# AUTHORS:      Bruce Van Aartsen
# AFFILIATION:  NASA GSFC / SSAI
# DATE:         Feb 1, 2017
#
# DESCRIPTION:
# This script will parse *.cnt output files from the MET GridStat script to 
# graph a time-series of RMSE, BIAS, and Pearson Correlation for given parameters
# 
# This script is meant to be executed by Jupyter Notebook: NUWRF_validation.ipynb
# It requires the following parameters to be predefined in the Notebook environment:
# fileDate - Python list of analysis periods in this format: YYYYMMDDHH
# inputRootDir - Directory that contains the MET PointStat output *.cnt files
# nuwrfRunDirs - Python list of subdirectories containing NU-WRF ensemble runs
# nuwrfRunLabels - Python list of labels corresponding to above ensemble runs
#                                   
#------------------------------------------------------------------------------

import os, sys
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

#beginYear = 2014
#beginMon = 04
#beginDate = 28
#beginHour = 00
#endYear = 2014
#endMon = 05
#endDate = 03
#endHour = 00
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
# Set up string variables for the input files
#------------------------------------------------------------------------------

startFilename = "grid_stat_APCP_24_"
midFilename = "0000L_"
endFilename = "0000V_cnt.txt"
#inputRootDir = "/discover/nobackup/projects/nu-wrf/members/bvanaart/SkillScoreRuns/MET_output/GridStat/round"
#exptDirList = ["0_default", "1_mp_physics_6","1_mp_physics_8","1_mp_physics_10","1_mp_physics_16","1_mp_physics_17","1_mp_physics_19","1_mp_physics_21","1_mp_physics_55"]
exptDirList = nuwrfRunDirs
#exptLabels = [" Def_Goddard4"," WSM6_graupel", " Thomson 2-m", " Morrison 2-m"," WDM6 2-m_gr"," NSSL 2-m"," NSSL 1-m"," NSSL-LFO 1-m"," Goddard 3ICE"]
exptLabels = nuwrfRunLabels
myColors = [colorConverter.to_rgba(c)
    for c in ('red','black','goldenrod','blue','darkseagreen','darkgreen','mediumpurple','indigo','darkturquoise')]
numExpts = len(exptLabels)

#Set output PDF file name, which will host several plots
eLabel = '-'.join(str(x) for x in exptLabels[0:2])
pdfName = 'Precip-RMSE-Bias_' + ''.join(eLabel.split()) + '.pdf'
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
#print dateIncr
allDates = np.array([startDate + datetime.timedelta(hours = int(hourIncr*i)) for i in range(numPeriods)])
#Double check:
for dateTime in allDates:
    print (dateTime.strftime("%Y-%m-%d-%HZ"))

#------------------------------------------------------------------------------
# Set up arrays to hold stats values
#------------------------------------------------------------------------------

biasPrecip = np.zeros((numExpts,numPeriods))
rmsePrecip = np.zeros((numExpts,numPeriods))
corrPrecip = np.zeros((numExpts,numPeriods))

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

      #------------------------------------------------------------------------
      # Parse *cnt file to extract desired stats. Most useful:
      # Column#  Name         Description
      # -------  ----         -----------
      # 9        FCST_VAR     Model variable (i.e., T2, Q2, U10, V10)
      # 15       INTERP_MTHD  Interpolation method (NEAREST, MEDIAN, DW_MEAN)
      # 22       TOTAL        Total number of matched pairs for comparison
      # 28       FSTDEV       Standard deviation of the forecasts
      # 38       OSTDEV       Standard deviation of the observations
      # 43       PR_CORR      Pearson correlation coefficient
      # 53       ME           Mean Error (AKA: Bias)
      # 66       MAE          Mean absolute error
      # 75       RMSE         Root mean squared error
      #------------------------------------------------------------------------

      file = open(pathName, 'r')
      # Use csv.reader to read each line and split it into columns: line[col#]
      for line in csv.reader(file, delimiter=" ", skipinitialspace=True):
         # Check column 20 for ALPHA=0.1 and col 14 for VX_MASK=WRF_Ens
         # where line index = (column# - 1)
         if (line[19] == "0.1") and (line[13] == "WRF_Ens"):
            
            # Get RMSE and Bias from lines that contain RAINNC0_24
            if (line[8] == "RAINNC0_24"):
               biasPrecip[exptNum,periodNum] = float(line[52])
               rmsePrecip[exptNum,periodNum] = float(line[74])
               corrPrecip[exptNum,periodNum] = float(line[42])
               
            
#print rmsePrecip

#------------------------------------------------------------------------------
# Start plotting the data
#------------------------------------------------------------------------------

#Create doc to host several plots (now defined near top)
#pp = PdfPages('Precip-RMSE_Microphys.pdf')

#Precip plots
#fig = plt.figure()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, rmsePrecip[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('24-hr Precipitation RMSE (mm)')
plt.title('24-hr Precip RMSE Time Series ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Precip_RMSE.png')

#fig = plt.figure()
plt.clf()

for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, biasPrecip[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

plt.ylabel('24-hr Precipitation Bias (mm)')
plt.title('24-hr Precip Bias Time Series ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Precip_Bias.png')


#Pearson correlation coefficient plot
#fig = plt.figure()
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, corrPrecip[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
plt.ylabel('24-hr Precipitation Correlation')
plt.title('Pearson Precip Correlation Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('Precip_Pearson.png')


#Close PDF doc
pp.close()
