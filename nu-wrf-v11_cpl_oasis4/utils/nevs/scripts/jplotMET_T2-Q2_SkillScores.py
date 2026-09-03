#!/usr/bin/env python
#------------------------------------------------------------------------------
# AUTHORS:      Bruce Van Aartsen
# AFFILIATION:  NASA GSFC / SSAI
# DATE:         May 27, 2016
#
# DESCRIPTION:
# This script will parse *.cts output files from the MET PointStat script to 
# graph a time-series of Skill Scores for given parameters
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

startFilename = "point_stat_"
midFilename = "0000L_"
endFilename = "0000V_cts.txt"
#inputRootDir = "/discover/nobackup/projects/nu-wrf/members/bvanaart/SkillScoreRuns/MET_output/PointStat/round"
#exptDirList = ["0_default","2_bl_pbl_1_sf_sfclay_11","2_bl_pbl_2_sf_sfclay_2","3_ra_lw_phys_4_sw_phys_4","3_ra_lw_phys_56_sw_phys_56","4_moist_adv_opt_2","4_moist_adv_opt_4"]
exptDirList = nuwrfRunDirs

#exptLabels = [" Def_Goddard4"," YSU PBL, MM5 M-O SL", " MYJ TKE PBL, M-O-J SL", " RRTMG"," 2014 Goddard rad"," Monotonic advec"," 5th-order WENO advec"]
exptLabels = nuwrfRunLabels
myColors = [colorConverter.to_rgba(c)
    for c in ('red','black','goldenrod','blue','darkseagreen','darkgreen','mediumpurple','indigo','darkturquoise')]
numExpts = len(exptLabels)

#Set output PDF file name, which will host several plots
eLabel = '-'.join(str(x) for x in exptLabels[0:2])
pdfName = 'T-Q-U-V_SkillScores_' + ''.join(eLabel.split()) + '.pdf'
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
print ('Time increment = ' + str(dateIncr))
allDates = np.array([startDate + datetime.timedelta(hours = int(hourIncr*i)) for i in range(numPeriods)])
#Double check:
for dateTime in allDates:
    print (dateTime.strftime("%Y-%m-%d-%HZ"))

#------------------------------------------------------------------------------
# Set up arrays to hold stats values
#------------------------------------------------------------------------------

GSST2 = np.zeros((numExpts,numPeriods))
HKT2 = np.zeros((numExpts,numPeriods))
HSST2 = np.zeros((numExpts,numPeriods))
GSSQ2 = np.zeros((numExpts,numPeriods))
HKQ2 = np.zeros((numExpts,numPeriods))
HSSQ2 = np.zeros((numExpts,numPeriods))
GSSU10 = np.zeros((numExpts,numPeriods))
HKU10 = np.zeros((numExpts,numPeriods))
HSSU10 = np.zeros((numExpts,numPeriods))
GSSV10 = np.zeros((numExpts,numPeriods))
HKV10 = np.zeros((numExpts,numPeriods))
HSSV10 = np.zeros((numExpts,numPeriods))

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
# 17       FCST_THRESH  Threshold value to test for correct forecast
# 22       TOTAL        Total number of matched pairs for comparison
# 66       GSS          Gilbert Skill Score
# 69       HK           Hanssen-Kuipers Discriminant (True skill statistic)
# 74       HSS          Heidke Skill Score
#------------------------------------------------------------------------------

      file = open(pathName, 'r')
      # Use csv.reader to read each line and split it into columns
      for line in csv.reader(file, delimiter=" ", skipinitialspace=True):
         # Check column 15 to see if INTERP_MTHD is DW_MEAN
         if (line[14] == "DW_MEAN"):
            
            # Get Skill Scores from lines where T2 >300 (use column# - 1)
            if (line[8] == "T2") and (line[16] == ">300.0"):
               GSST2[exptNum,periodNum] = float(line[65])
               HKT2[exptNum,periodNum] = float(line[68])
               HSST2[exptNum,periodNum] = float(line[73])
               
            # Get Skill Scores from lines where Q2 >0.016 (use column# - 1)
            if (line[8] == "Q2") and (line[16] == ">0.016"):
               GSSQ2[exptNum,periodNum] = float(line[65])
               HKQ2[exptNum,periodNum] = float(line[68])
               HSSQ2[exptNum,periodNum] = float(line[73])
            
            # Get Skill Scores from lines where U10 >5 m/s (use column# - 1)
            if (line[8] == "U10") and (line[16] == ">5.0"):
               #print line[65],line[68],line[73]
               GSSU10[exptNum,periodNum] = float(line[65]) if line[65] != "NA" else  0.0
               HKU10[exptNum,periodNum] = float(line[68]) if line[68] != "NA" else  0.0
               HSSU10[exptNum,periodNum] = float(line[73]) if line[73] != "NA" else  0.0
               
            # Get Skill Scores from lines where V10 >5 m/s (use column# - 1)
            if (line[8] == "V10") and (line[16] == ">5.0"):
               GSSV10[exptNum,periodNum] = float(line[65]) if line[65] != "NA" else  0.0
               HKV10[exptNum,periodNum] = float(line[68]) if line[68] != "NA" else  0.0
               HSSV10[exptNum,periodNum] = float(line[73]) if line[73] != "NA" else  0.0
            
#print GSST2
#print HKT2
#      print "rmseT2 = %f" %rmseT2[0,0]
#      print "corrT2 = %f" %corrT2[0,0]
#      print "biasQ2 = %f" %biasQ2[0,0]
#      print "rmseQ2 = %f" %rmseQ2[0,0]
#      print "corrQ2 = %f" %corrQ2[0,0]

#------------------------------------------------------------------------------
# Start plotting the data
#------------------------------------------------------------------------------

#Create doc to host several plots (now defined near top)
#pp = PdfPages('T-Q-U-V_SkillScores.pdf')

#-----------------------------
### Temperature (T2) plots ###
#-----------------------------
# T2 Gilbert Skill
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, GSST2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'Gilbert Skill (T2 > 27\u00B0C)')
plt.title(u'Gilbert Skill Score - Temp > 27\u00B0C ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Temp-2m_GilbertSkillScore.png')

# T2 HK True Skill
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HKT2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'HK Discriminant (T2 > 27\u00B0C)')
plt.title(u'Hanssen-Kuipers (True Skill) - Temp > 27\u00B0C')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Temp-2m_H-KDiscriminant.png')

# T2 Heidke Skill
plt.clf()

for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HSST2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'Heidke Skill Score (T2 > 27\u00B0C)')
plt.title(u'Heidke Skill Score - Temp > 27\u00B0C')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('Temp-2m_HeidkeSkillScore.png')

#-----------------------------
#Specific Humidity (Q2) plots
#-----------------------------
# Q2 Gilbert Skill
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, GSSQ2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'Gilbert Skill (Q2 > 0.016)')
plt.title(u'Gilbert Skill Score - Spec Humidity > 0.016 ')
plt.legend(loc='lower left', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Spec-Humidity_GilbertSkillScore.png')

# Q2 HK True Skill
plt.clf()

for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HKQ2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'HK Discriminant (Q2 > 0.016)')
plt.title(u'Hanssen-Kuipers (True Skill) - Spec Humidity > 0.016')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('Spec-Humidity_H-KDiscriminant.png')

# Q2 Heidke Skill
plt.clf()

for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HSSQ2[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'Heidke Skill Score (Q2 > 0.016)')
plt.title(u'Heidke Skill Score - Spec Humidity > 0.016 ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()

#-----------------------------
### U-Wind (U10) plots ###
#-----------------------------
# U10 Gilbert Skill
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, GSSU10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'Gilbert Skill (U10 > 5 m/s)')
plt.title(u'Gilbert Skill Score - U-Wind > 5 m/s ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('U-Wind_GilbertSkillScore.png')

# U10 HK True Skill
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HKU10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'HK Discriminant (U10 > 5 m/s)')
plt.title(u'Hanssen-Kuipers (True Skill) - U-Wind > 5 m/s ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('U-Wind_H-KDiscriminant.png')

# U10 Heidke Skill
plt.clf()

for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HSSU10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'Heidke Skill Score (U10 > 5 m/s)')
plt.title(u'Heidke Skill Score - U-Wind > 5 m/s ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('U-Wind_HeidkeSkillScore.png')

#-----------------------------
# V-Wind (V10) plots
#-----------------------------
# V10 Gilbert Skill
plt.clf()
for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, GSSV10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'Gilbert Skill (V10 > 5 m/s)')
plt.title(u'Gilbert Skill Score - V-Wind > 5 m/s ')
plt.legend(loc='lower left', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('V-Wind_GilbertSkillScore.png')

# V10 HK True Skill
plt.clf()

for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HKV10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()

#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'HK Discriminant (V10 > 5 m/s)')
plt.title(u'Hanssen-Kuipers (True Skill) - V-Wind > 5 m/s ')
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show()
pp.savefig()
pylab.savefig('V-Wind_H-KDiscriminant.png')

# V10 Heidke Skill
plt.clf()

for exptNum in range(numExpts):
    plt.gca().xaxis.set_major_locator(DayLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter('%m/%d-%HZ'))
    plt.plot_date(allDates, HSSV10[exptNum,:], 'o-', label=exptLabels[exptNum], color=myColors[exptNum], linewidth=2)
    plt.gcf().autofmt_xdate()
    
#Use this to get min/max of plot and then reset y-range to include 0    
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,y2))

plt.ylabel(u'Heidke Skill Score (V10 > 5 m/s)')
plt.title(u'Heidke Skill Score - V-Wind > 5 m/s ')    
plt.legend(loc='best', shadow=True, fontsize='x-small', fancybox=True)
#plt.show() 
pp.savefig()
pylab.savefig('V-Wind_HeidkeSkillScore.png')

#Close PDF doc
pp.close()
plt.close()
