#!/usr/bin/env python
#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Sciences and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# SCRIPT: 
# c1440_NR.geos2wrf.py
#
# AUTHOR:
# Eric Kemp, NASA CISTO/SSAI
# 
# DESCRIPTION:
# Script for processing GEOS-5 7-km Ganymed Nature Run
#
#------------------------------------------------------------------------------

# Standard modules
import ConfigParser
import datetime
import os
import subprocess
import sys

# Other modules
import geos2wrfUtil

# Process command line arguments
if len(sys.argv) != 3:
    text  = "ERROR with command line arguments!"
    text += "/nUsage: %s settingsCfgFile variablesCfgFile" %(sys.argv[0])
    raise SystemExit, text

settingsCfgFile = sys.argv[1]
if not os.path.exists(settingsCfgFile):
    raise SystemExit,"ERROR, %s does not exist!" %(settingsCfgFile)
variablesCfgFile = sys.argv[2]
if not os.path.exists(variablesCfgFile):
    raise SystemExit,"ERROR, %s does not exist!" %(variablesCfgFile)

# Flag indicating if const file has been processed yet.
constDone=False

# Read settings config file. Data are save in global module variables.
g2w = geos2wrfUtil.Geos2Wrf(settingsCfgFile,variablesCfgFile)

# For each time, process each collection and then derive missing variables.
currentDateTime = g2w.startDateTime
while currentDateTime <= g2w.endDateTime:
    g2w.runGEOS2WPS(currentDateTime)

    # All collections processed at this point. Now derive missing variables.
    if g2w.needTerrain:
        if not g2w.constDone:
            print "Deriving terrain height..."
            g2w.runCreateSOILHGT(currentDateTime)

    if g2w.needLandSea:
        if not g2w.constDone:
            print "Deriving land sea mask..."
            g2w.runCreateLANDSEA(currentDateTime)

    if g2w.needRH:
        print "Deriving relative humidity..."
        g2w.runCreateRH(currentDateTime)

    # Consolidate and clean up
    g2w.consolidate(currentDateTime)

    # Move on to next time
    g2w.setConstDone(True)
    currentDateTime += g2w.timeDelta


