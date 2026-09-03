#!/usr/bin/env python

#------------------------------------------------------------------------------
#
# Define class to automate most invocations of GEOS2WRF software.
# WARNING: Customized mostly for GEOS-5 7-km Nature Run.
#
#------------------------------------------------------------------------------

class Geos2Wrf:
    """Class for automating GEOS2WRF package."""

    #--------------------------------------------------------------------------
    def __init__(self,settingsCfgFile,variablesCfgFile,
                 constDone=False):
        """Create object"""
        self.readSettingsCfgFile(settingsCfgFile)
        self.readVariablesCfgFile(variablesCfgFile)

        self.constDone = constDone

    #--------------------------------------------------------------------------
    def readSettingsCfgFile(self,settingsCfgFile):
        """Process config file with settings."""

        # Standard modules
        import ConfigParser

        # Build object to read the config file.
        config = ConfigParser.RawConfigParser()
        config.read(settingsCfgFile)

        self._readTimesSection(config)
        self._readPathsSection(config)
        self._readDeriveSection(config)
        self._readSubsetSection(config)

    #--------------------------------------------------------------------------
    def _readTimesSection(self,config):

        # Standard modules
        import datetime

        startTimeString = config.get("TIMES","startDateTime")
        minute = 0
        second = 0
        if len(startTimeString) == 14:
            second = int(startTimeString[12:14])
            minute = int(startTimeString[10:12])
        elif len(startTimeString) == 12:
            minute = int(startTimeString[10:12])
        elif len(startTimeString) != 10:
            raise SystemExit, "ERROR parsing startDateTime %s" \
                %(startTimeString)
        self.startDateTime = datetime.datetime(year=int(startTimeString[0:4]), 
                                               month=int(startTimeString[4:6]),
                                               day=int(startTimeString[6:8]),
                                               hour=int(startTimeString[8:10]),
                                               minute=minute,second=second)

        endTimeString = config.get("TIMES","endDateTime")
        if len(endTimeString) == 14:
            second = int(endTimeString[12:14])
            minute = int(endTimeString[10:12])
        elif len(endTimeString) == 12:
            minute = int(endTimeString[10:12])
        elif len(endTimeString) != 10:
            raise SystemExit, "ERROR parsing endDateTime %s" \
                %(endTimeString)
        self.endDateTime = datetime.datetime(year=int(endTimeString[0:4]), 
                                             month=int(endTimeString[4:6]),
                                             day=int(endTimeString[6:8]),
                                             hour=int(endTimeString[8:10]),
                                             minute=minute,second=second)

        deltaSeconds = config.getint("TIMES","deltaSeconds")
        self.timeDelta = datetime.timedelta(seconds=deltaSeconds)
        tavgOffsetSeconds = config.getint("TIMES","tavgOffsetSeconds")
        self.tavgTimeDelta = datetime.timedelta(seconds=tavgOffsetSeconds)

    #--------------------------------------------------------------------------
    def _readPathsSection(self,config):
        self.topGeosDataDir = config.get("PATHS","topGeosDataDir")
        self.geos2wps       = config.get("PATHS","geos2wps")
        self.createSOILHGT  = config.get("PATHS","createSOILHGT")
        self.createLANDSEA  = config.get("PATHS","createLANDSEA")
        self.createRH       = config.get("PATHS","createRH")

    #--------------------------------------------------------------------------
    def _readDeriveSection(self,config):
        self.needTerrain    = config.getboolean("DERIVE","needTerrain")
        self.needLandSea    = config.getboolean("DERIVE","needLandSea")
        self.needRH         = config.getboolean("DERIVE","needRH")

    #--------------------------------------------------------------------------
    def _readSubsetSection(self,config):
        self.spatialSubset    = config.getboolean("SUBSET","spatialSubset")
        self.iLonMin          = config.getint("SUBSET","iLonMin")
        self.iLonMax          = config.getint("SUBSET","iLonMax")
        self.jLatMin          = config.getint("SUBSET","jLatMin")
        self.jLatMax          = config.getint("SUBSET","jLatMax")
        self.kVertMin         = config.getint("SUBSET","kVertMin")
        self.kVertMax         = config.getint("SUBSET","kVertMax")

    #--------------------------------------------------------------------------
    def readVariablesCfgFile(self,variablesCfgFile):
        """Read config file with GEOS variables and store entries as data 
        attributes."""

        # Standard modules
        import ConfigParser

        # Build object to read the config file.
        config = ConfigParser.RawConfigParser()
        config.read(variablesCfgFile)

        # Store the entries as dictionary of dictionaries
        self.variables = {}
        sections=config.sections()
        for section in sections:
            if section not in self.variables.keys():
                self.variables[section] = {}
                for setting in ["collection","wpsName","units","description"]:
                    self.variables[section][setting] = \
                        config.get(section,setting)
                for setting in ["rank","levelType"]:
                    self.variables[section][setting] = \
                        config.getint(section,setting)

        # Also convenient to group variables by collection
        self.variablesByCollection = {}
        for section in sections:
            collection=config.get(section,"collection")
            if collection not in self.variablesByCollection.keys():
                self.variablesByCollection[collection] = []
            self.variablesByCollection[collection].append(section)

    #--------------------------------------------------------------------------
    def runGEOS2WPS(self,currentDateTime):
        """Method for running GEOS2WPS"""

        # Standard modules
        import os
        import subprocess

        cyear = currentDateTime.year
        cmonth = currentDateTime.month
        cday = currentDateTime.day
        chour = currentDateTime.hour
        cminute = currentDateTime.minute
        csecond = currentDateTime.second

        # Make sure time stamp reflects instantaneous time. Follow logic
        # in METGRID to determine whether minutes and seconds should be
        # included in the timeStamp.
        timeStamp = self._calcTimeStamp(currentDateTime)
        
        # Loop through and process each collection
        collections = self.variablesByCollection.keys()
        collections.sort()
        for collection in collections:
            numberOfVariables = len(self.variablesByCollection[collection])

            if "const" in collection:
                type = "const"
                # We only need to process the first const file
                if self.constDone:
                    print "Skipping ", collection
                    continue
            elif "inst" in collection:
                type = "inst"
            elif "tavg" in collection:
                type = "tavg"
            else:
                raise SystemExit, \
                    "ERROR, unknown collection type %s" %(collection)

            # Must fall back on earlier time if using tavg data
            if type == "tavg":
                tavgCurrentDateTime = currentDateTime - self.tavgTimeDelta
                cyear   = tavgCurrentDateTime.year
                cmonth  = tavgCurrentDateTime.month
                cday    = tavgCurrentDateTime.day
                chour   = tavgCurrentDateTime.hour
                cminute = tavgCurrentDateTime.minute

            # FIXME: Reorganize logic to make this easier to accomodate
            # other directory and file name conventions.
            path  = "%s"  %(self.topGeosDataDir)
            path += "/%s" %(type)
            path += "/%s" %(collection)
            path += "/Y%4.4d" %(cyear)
            path += "/M%2.2d" %(cmonth)
            path += "/D%2.2d" %(cday)
            if type == "const":
                path += "/c1440_NR.%s.%4.4d%2.2d%2.2d.nc4" \
                    %(collection,cyear,cmonth,cday)
            else:
                path += "/c1440_NR.%s.%4.4d%2.2d%2.2d_%2.2d%2.2dz.nc4" \
                    %(collection,cyear,cmonth,cday,chour,cminute)
            print "path = %s" %(path)

            # Create temporary symbolic link
            if os.path.lexists('GEOSFILE'):
                os.unlink('GEOSFILE')
            print 'Creating symbolic link: GEOSFILE -> %s' %(path)
            os.symlink(path,'GEOSFILE')
            if not os.path.exists('GEOSFILE'):
                raise SystemExit, "ERROR, cannot link to %s" %(path)

            # Create temporary namelist file for GEOS2WPS
            fd = open("namelist.geos2wps","w")
            fd.write("""
&files
   geosFileFormat = 2,
   geosFileName = 'GEOSFILE',
   outputDirectory = './',
/
""")

            fd.write("""
&coordinates
   longitudeName='lon',
   latitudeName='lat',
   hasVerticalDimension=.true.,
   verticalName='lev',
/
""")

            fd.write("""
&forecast
   numberOfTimes=1,
   validTimes(1)='%s'
   timeIndices(1)=1,
   forecastHours(1)=0,
/
""" \
                         %(timeStamp))

            fd.write("""
&variables
   numberOfVariables=%d,
""" \
                         %(numberOfVariables))

            for i in range(1,numberOfVariables+1):
                section = self.variablesByCollection[collection][i-1]
                fd.write("""
   variableRanks(%d) = %d,
   variableLevelTypes(%d) = %d,
   variableNamesIn(%d) = '%s',
   variableNamesOut(%d) = '%s',
   variableUnits(%d) = '%s',
   variableDescriptions(%d) = '%s',
""" \
                             %(i, int(self.variables[section]["rank"]),
                               i, int(self.variables[section]["levelType"]),
                               i, section,
                               i, self.variables[section]["wpsName"],
                               i, self.variables[section]["units"],
                               i, self.variables[section]["description"]))

            
            fd.write("""
/

""")

            # Only write this block if subset is needed. GEOS2WPS 
            # processes entire GEOS grid by default.
            if self.spatialSubset:                
                fd.write("""

&subsetData
   subset=.true.,
   iLonMin=%d,
   iLonMax=%d,
   jLatMin=%d,
   jLatMax=%d,
   kVertMin=%d,
   kVertMax=%d,
/
""" \
                             %(self.iLonMin,self.iLonMax,
                               self.jLatMin,self.jLatMax,
                               self.kVertMin,self.kVertMax))
            else:
                fd.write("""

&subsetData
   subset=.false.,
/
""")
            # Make sure file is completely written
            fd.close()

            # Run GEOS2WPS
            print "Processing", collection
            cmd = ['%s' %(self.geos2wps)]
            try:
                subprocess.check_call(cmd)
            except:
                raise SystemExit,'ERROR calling geos2wps!'

            # Clean up
            os.unlink('GEOSFILE')

            # Go to next collection
            continue

    #--------------------------------------------------------------------------
    def runCreateSOILHGT(self,currentDateTime):
        """Method for running createSOILHGT"""

        # Standard modules
        import os
        import shutil
        import subprocess

        timeStamp = self._calcTimeStamp(currentDateTime)
        includeMinute,includeSecond = self._setIncludeMinuteSecond()

        # Create temporary WPS file
        src = "PHIS_GROUND_LEVEL:%s" %(timeStamp)
        dst = "GEOS_TMP:%s" %(timeStamp)
        shutil.copyfile(src,dst)

        # Create temporary namelist file for createSOILHGT
        fd = open("namelist.createSOILHGT","w")
        fd.write("""
&input
   directory='./',
   prefix='GEOS_TMP',
   year=%4.4d,
   month=%2.2d,
   day=%2.2d,
   hour=%2.2d,
   minute=%2.2d,
   second=%2.2d,
   includeMinute=%s,
   includeSecond=%s,
   surfaceGeopotentialName='PHIS',
/

""" %(currentDateTime.year,
      currentDateTime.month,
      currentDateTime.day,
      currentDateTime.hour,
      currentDateTime.minute,
      currentDateTime.second,
      includeMinute,includeSecond))
        fd.close()

        # Run createSOILHGT
        cmd = ['%s' %(self.createSOILHGT)]
        try:
            subprocess.check_call(cmd)
        except:
            raise SystemExit,'ERROR calling createSOILHGT!'

        # Clean up
        tmpfile = "GEOS_TMP:%s" %(timeStamp)
        os.unlink(tmpfile)

        # Terrain height is time invariant
        src = "SOILHGT_GROUND_LEVEL:%s" %(timeStamp)
        dst = "SOILHGT_GROUND_LEVEL"
        shutil.move(src,dst)

    #--------------------------------------------------------------------------
    def runCreateLANDSEA(self,currentDateTime):
        """Method for running createLANDSEA"""

        # Standard modules
        import os
        import shutil
        import subprocess

        timeStamp = self._calcTimeStamp(currentDateTime)
        includeMinute,includeSecond = self._setIncludeMinuteSecond()

        # Cat the FRLAKE and FROCEAN data to a new temporary file
        frlake  = "FRLAKE_GROUND_LEVEL:%s" %(timeStamp)
        frocean = "FROCEAN_GROUND_LEVEL:%s" %(timeStamp)
        geos_tmp = "GEOS_TMP:%s" %(timeStamp)
        if os.path.exists(geos_tmp):
            os.unlink(geos_tmp)
        # Having problems with shell redirection with subprocess. So 
        # we'll fall back on "obsolescent" approach.
        os.system("cat %s %s > %s" %(frlake,frocean,geos_tmp))

        # Create temporary namelist file for createLANDSEA
        fd = open("namelist.createLANDSEA","w")
        fd.write("""
&input
   directory='./',
   prefix='GEOS_TMP',
   year=%4.4d,
   month=%2.2d,
   day=%2.2d,
   hour=%2.2d,
   minute=%2.2d,
   second=%2.2d,
   includeMinute=%s,
   includeSecond=%s,
   lakeFractionName='FRLAKE',
   oceanFractionName='FROCEAN',
/

""" %(currentDateTime.year,
      currentDateTime.month,
      currentDateTime.day,
      currentDateTime.hour,
      currentDateTime.minute,
      currentDateTime.second,
      includeMinute,includeSecond))
        fd.close()

        # Run createLANDSEA
        cmd = ['%s' %(self.createLANDSEA)]
        try:
            subprocess.check_call(cmd)
        except:
            raise SystemExit,'ERROR calling createLANDSEA!'

        # Clean up
        tmpfile = "GEOS_TMP:%s" %(timeStamp)
        os.unlink(geos_tmp)

        # Terrain height is time invariant
        src = "LANDSEA_GROUND_LEVEL:%s" %(timeStamp)
        dst = "LANDSEA_GROUND_LEVEL"
        shutil.move(src,dst)

    #--------------------------------------------------------------------------
    def runCreateRH(self,currentDateTime):
        """Runs createRH"""

        # Standard modules
        import os
        import subprocess

        timeStamp = self._calcTimeStamp(currentDateTime)
        includeMinute,includeSecond = self._setIncludeMinuteSecond()

        # Cat the input variables to a new temporary file
        psfc  = "PSFC_GROUND_LEVEL:%s" %(timeStamp)
        pressure = "PRESSURE_MODEL_LEVEL:%s" %(timeStamp)
        tt_model = "TT_MODEL_LEVEL:%s" %(timeStamp)
        tt_2m = "TT_2M_ABOVE_GROUND_LEVEL:%s" %(timeStamp)
        spechumd_model = "SPECHUMD_MODEL_LEVEL:%s" %(timeStamp)
        spechumd_2m = "SPECHUMD_2M_ABOVE_GROUND_LEVEL:%s" %(timeStamp)
        geos_tmp = "GEOS_TMP:%s" %(timeStamp)
        if os.path.exists(geos_tmp):
            os.unlink(geos_tmp)
            # Having problems with shell redirection with subprocess. So we'll
            # fall back on "obsolescent" approach.
        os.system("cat %s %s %s %s %s %s > %s" \
                  %(psfc,pressure,tt_model,tt_2m,spechumd_model,spechumd_2m,
                    geos_tmp))

        # Create temporary namelist file for createLANDSEA
        fd = open("namelist.createRH","w")
        fd.write("""
&input
   directory='./',
   prefix='GEOS_TMP',
   year=%4.4d,
   month=%2.2d,
   day=%2.2d,
   hour=%2.2d,
   minute=%2.2d,
   second=%2.2d,
   includeMinute=%s,
   includeSecond=%s,
   processSurfacePressure=.true.,
   surfacePressureName='PSFC',
   pressureName='PRESSURE',
   temperatureName='TT',
   specificHumidityName='SPECHUMD',
/

""" %(currentDateTime.year,
      currentDateTime.month,
      currentDateTime.day,
      currentDateTime.hour,
      currentDateTime.minute,
      currentDateTime.second,
      includeMinute,includeSecond))
        fd.close()

        # Run createRH
        cmd = ['%s' %(self.createRH)]
        try:
            subprocess.check_call(cmd)
        except:
            raise SystemExit,'ERROR calling createRH!'

        # Clean up
        tmpfile = "GEOS_TMP:%s" %(timeStamp)
        os.unlink(geos_tmp)

    #--------------------------------------------------------------------------
    def consolidate(self,currentDateTime):
        """Consolidate GEOS data for current date and time to single file."""

        # Standard modules
        import os
        import shutil

        timeStamp = self._calcTimeStamp(currentDateTime)

        # Cat all data for WPS into single file.
        soilhgt   = "SOILHGT_GROUND_LEVEL"
        landsea   = "LANDSEA_GROUND_LEVEL"
        hgt       = "HGT_MODEL_LEVEL:%s"           %(timeStamp)
        pmsl      = "PMSL_MEAN_SEA_LEVEL:%s"       %(timeStamp)
        pressure  = "PRESSURE_MODEL_LEVEL:%s"      %(timeStamp)
        psfc      = "PSFC_GROUND_LEVEL:%s"         %(timeStamp)
        rh        = "RH_MODEL_LEVEL:%s"            %(timeStamp)
        rh_2m     = "RH_2M_ABOVE_GROUND_LEVEL:%s"  %(timeStamp)
        seaice    = "SEAICE_GROUND_LEVEL:%s"       %(timeStamp)
        skintemp  = "SKINTEMP_GROUND_LEVEL:%s"     %(timeStamp)
        tt_2m     = "TT_2M_ABOVE_GROUND_LEVEL:%s"  %(timeStamp)
        tt        = "TT_MODEL_LEVEL:%s"            %(timeStamp)
        uu_10m    = "UU_10M_ABOVE_GROUND_LEVEL:%s" %(timeStamp)
        uu        = "UU_MODEL_LEVEL:%s"            %(timeStamp)
        vv_10m    = "VV_10M_ABOVE_GROUND_LEVEL:%s" %(timeStamp)
        vv        = "VV_MODEL_LEVEL:%s"            %(timeStamp)
        geos_save = "GEOS:%s.save"                 %(timeStamp)
        geos      = "GEOS:%s"                      %(timeStamp)
        if os.path.exists(geos_save):
            os.unlink(geos_save)
        # Having problems with shell redirection with subprocess. So we'll
        # fall back on "obsolescent" approach.
        cmd = "cat " + 16*"%s " + "> %s"
        os.system(cmd \
                  %(landsea,pressure,psfc,pmsl,tt,tt_2m,skintemp,seaice,
                    rh,rh_2m,uu_10m,uu,vv_10m,vv,hgt,soilhgt,
                    geos_save))
    
        # Save some space
        print "Purging files..."
        cmd = "rm *:%s" %(timeStamp)
        print cmd
        os.system(cmd)

        # Now rename GEOS file
        src=geos_save
        dst=geos
        print "Renaming %s to %s" %(src,dst)
        shutil.move(src,dst)

    #--------------------------------------------------------------------------
    def setConstDone(self,constDone):
        """Setter for constDone"""
        self.constDone=constDone
    #--------------------------------------------------------------------------
    def _calcTimeStamp(self,currentDateTime):
        """Common internal method for making time stamp"""

        cyear = currentDateTime.year
        cmonth = currentDateTime.month
        cday = currentDateTime.day
        chour = currentDateTime.hour
        cminute = currentDateTime.minute
        csecond = currentDateTime.second

        includeMinute,includeSecond = self._setIncludeMinuteSecond()

        if includeMinute == ".true." and includeSecond == ".true.":
            timeStamp = "%4.4d-%2.2d-%2.2d_%2.2d:%2.2d:%2.2d" \
                %(cyear,cmonth,cday,chour,cminute,csecond)
        elif includeMinute == ".true.":
            timeStamp = "%4.4d-%2.2d-%2.2d_%2.2d:%2.2d" \
                %(cyear,cmonth,cday,chour,cminute)
        else:
            timeStamp = "%4.4d-%2.2d-%2.2d_%2.2d" %(cyear,cmonth,cday,chour)
        
        return timeStamp

    #--------------------------------------------------------------------------
    def _setIncludeMinuteSecond(self):
        """Returns text on whether WPS file names should include minutes
        and seconds"""

        if self.timeDelta.seconds % 3600 == 0:
            includeMinute=".false."
            includeSecond=".false."
        elif self.timeDelta.seconds % 60 == 0:
            includeMinute=".true."
            includeSecond=".false."
        else:            
            includeMinute=".true."
            includeSecond=".true."

        return includeMinute,includeSecond
    #--------------------------------------------------------------------------
