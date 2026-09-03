#!/usr/bin/env python
#------------------------------------------------------------------------------
# AUTHORS:      Bruce Van Aartsen
# AFFILIATION:  NASA GSFC / SSAI
# DATE:         Feb 1, 2017
#
# DESCRIPTION:
# This script will parse *cnt.txt output files from the MET PointStat script to 
# graph a time-series of RMSE, BIAS, and Pearson Correlation for given parameters
# 
# This script is meant to be executed by Jupyter Notebook: NUWRF_validation.ipynb
# It requires the following parameters to be predefined in the Notebook environment:
# fileDate - Python list of analysis periods in this format: YYYYMMDDHH
# inputRootDir - Directory that contains the MET PointStat output *cnt.txt files
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

startFilename = "point_stat_"
midFilename = "0000L_"
endFilename = "0000V_cnt.txt"
#inputRootDir = "/discover/nobackup/projects/nu-wrf/members/bvanaart/SkillScoreRuns/MET_output/PointStat/round" 
#exptDirList = ["0_default", "1_mp_physics_6","1_mp_physics_8","1_mp_physics_10","1_mp_physics_16","1_mp_physics_17","1_mp_physics_19","1_mp_physics_21","1_mp_physics_55"]
exptDirList = nuwrfRunDirs
print (exptDirList)
#exptLabels = [" Def_Goddard4"," WSM6_graupel", " Thomson 2-m", " Morrison 2-m"," WDM6 2-m_gr"," NSSL 2-m"," NSSL 1-m"," NSSL-LFO 1-m"," Goddard 3ICE"]
exptLabels = nuwrfRunLabels
myColors = [colorConverter.to_rgba(c)
    for c in ('red','black','goldenrod','blue','darkseagreen','darkgreen','mediumpurple','indigo','darkturquoise')]
numExpts = len(exptLabels)
#print numExpts

#Set output PDF file name, which will host several plots
eLabel = '-'.join(str(x) for x in exptLabels[0:2])
pdfName = 'T-Q-U-V_RMSE-Bias_' + ''.join(eLabel.split()) + '.pdf'
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
#numPeriods = len(fileDate) -1
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

biasT2 = np.zeros((numExpts,numPeriods))
rmseT2 = np.zeros((numExpts,numPeriods))
corrT2 = np.zeros((numExpts,numPeriods))
biasQ2 = np.zeros((numExpts,numPeriods))
rmseQ2 = np.zeros((numExpts,numPeriods))
corrQ2 = np.zeros((numExpts,numPeriods))
biasU10 = np.zeros((numExpts,numPeriods))
rmseU10 = np.zeros((numExpts,numPeriods))
corrU10 = np.zeros((numExpts,numPeriods))
biasV10 = np.zeros((numExpts,numPeriods))
rmseV10 = np.zeros((numExpts,numPeriods))
corrV10 = np.zeros((numExpts,numPeriods))

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
#------------------------------------------------------------------------------

      file = open(pathName, 'r')
      # Use csv.reader to read each line and split it into columns
      for line in csv.reader(file, delimiter=" ", skipinitialspace=True):
         # Check column 15 to see if INTERP_MTHD is DW_MEAN
         if (line[14] == "DW_MEAN"):
            
            # Get RMSE and Bias from lines that contain T2 (use column# - 1)
            if (line[8] == "T2"):
               biasT2[exptNum,periodNum] = float(line[52])
               rmseT2[exptNum,periodNum] = float(line[74])
               corrT2[exptNum,periodNum] = float(line[42])
               
            # Get RMSE and Bias from lines that contain Q2 (use column# - 1)
            if (line[8] == "Q2"):
               biasQ2[exptNum,periodNum] = float(line[52])
               rmseQ2[exptNum,periodNum] = float(line[74])
               corrQ2[exptNum,periodNum] = float(line[42])
            
            # Get RMSE and Bias from lines that contain U10 (use column# - 1)
            if (line[8] == "U10"):
               biasU10[exptNum,periodNum] = float(line[52])
               rmseU10[exptNum,periodNum] = float(line[74])
               corrU10[exptNum,periodNum] = float(line[42])
               
            # Get RMSE and Bias from lines that contain V10 (use column# - 1)
            if (line[8] == "V10"):
               biasV10[exptNum,periodNum] = float(line[52])
               rmseV10[exptNum,periodNum] = float(line[74])
               corrV10[exptNum,periodNum] = float(line[42])
            
#print biasT2
#print biasQ2
#      print "rmseT2 = %f" %rmseT2[0,0]
#      print "corrT2 = %f" %corrT2[0,0]
#      print "biasQ2 = %f" %biasQ2[0,0]
#      print "rmseQ2 = %f" %rmseQ2[0,0]
#      print "corrQ2 = %f" %corrQ2[0,0]
#print "rmseU10 = ", rmseU10
#------------------------------------------------------------------------------
# Start plotting the data
#------------------------------------------------------------------------------

#Create doc to host several plots (now defined near top)
#pp = PdfPages('T-Q-U-V_RMSE_Bias.pdf')


#-----------------------------
### Temperature (T2) plots ###
#-----------------------------
# T2 RMSE
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, rmseT2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'T2 RMSE (\u00B0C) @2M')
plt.title(u'Temperature RMSE Time Series ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Temp-2m_RMSE.png')

# T2 Bias
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, biasT2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

plt.ylabel(u'T2 Bias (\u00B0C) @2M')
plt.title(u'Temperature Bias Time Series ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Temp-2m_Bias.png')

# T2 Pearson Correlation
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, corrT2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
plt.ylabel('T2 Correlation @2M')
plt.title('Pearson T2 Correlation Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('Temp-2m_Pearson.png')

#-----------------------------
#Specific Humidity (Q2) plots
#-----------------------------
# Q2 Bias plot
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, rmseQ2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('Q2 (kg water vapor/kg dry air) @2M')
plt.title('Specific Humidity RMSE  Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('Spec-Humidity_RMSE.png')

# Q2 Bias plot
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, biasQ2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
plt.ylabel('Q2 (kg water vapor/kg dry air) @2M')
plt.title('Specific Humidity Bias  Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('Spec-Humidity_Bias.png')

# Q2 Pearson correlation coefficient plot
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, corrQ2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
plt.ylabel('Q2 Correlation @2M')
plt.title('Pearson Q2 Correlation Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('Spec-Humidity_Pearson.png')

#-----------------------------
# U-Wind (U10) plots
#-----------------------------
# U10 Bias plot
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, rmseU10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('U-Wind RMSE (m/s) @10M')
plt.title('U-Wind RMSE  Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('U-Wind_RMSE.png')

# U10 Bias plot
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, biasU10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
plt.ylabel('U-Wind Bias(m/s) @10M')
plt.title('U-Wind Bias  Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('U-Wind_Bias.png')

# U10 Pearson correlation coefficient plot
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, corrU10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
plt.ylabel('U-Wind Correlation @10M')
plt.title('Pearson U-Wind Correlation Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('U-Wind_Pearson.png')

#-----------------------------
# V-Wind (V10) plots
#-----------------------------
# V10 Bias plot
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, rmseV10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel('V-Wind RMSE (m/s) @10M')
plt.title('V-Wind RMSE  Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('V-Wind_RMSE.png')

# V10 Bias plot
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, biasV10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
plt.ylabel('V-Wind Bias(m/s) @10M')
plt.title('V-Wind Bias  Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('V-Wind_Bias.png')

# V10 Pearson correlation coefficient plot
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, corrV10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
plt.ylabel('V-Wind Correlation @10M')
plt.title('Pearson V-Wind Correlation Time Series ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('V-Wind_Pearson.png')

#Close PDF doc
pp.close()
plt.close()
