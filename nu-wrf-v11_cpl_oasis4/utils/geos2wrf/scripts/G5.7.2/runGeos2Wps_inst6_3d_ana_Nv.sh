#!/bin/sh 
# Sample script for processing GEOS inst6_3d_ana_Nv file with geos2wps.

startTimeString=$1
endTimeString=$2
deltaHour=$3
cfgFile=$4

source $cfgFile

syear=`echo $startTimeString | cut -c1-4`
smonth=`echo $startTimeString | cut -c5-6`
sday=`echo $startTimeString | cut -c7-8`
shour=`echo $startTimeString | cut -c9-10`

sdate="${syear}${smonth}${sday}"

currentTimeString=${startTimeString}
while [ "$currentTimeString" -le "$endTimeString" ] ; do

    cyear=`echo $currentTimeString | cut -c1-4`
    cmonth=`echo $currentTimeString | cut -c5-6`
    cday=`echo $currentTimeString | cut -c7-8`
    chour=`echo $currentTimeString | cut -c9-10`
    timeStamp="${cyear}-${cmonth}-${cday}_${chour}"

    cdate="${cyear}${cmonth}${cday}"

    timeStamp1=$timeStamp

    currentTimeString="${cyear}${cmonth}${cday}${chour}"
    forecastHour=`${calcForecastHour} $startTimeString $currentTimeString`

    cat > namelist.geos2wps<<EOF
&files
  geosFileFormat=2,
  geosFileName='e572p5_fp.inst3_3d_asm_Nv.${sdate}_${shour}z+${cdate}_${chour}00z.nc4'
  outputDirectory='./',
/

&coordinates
 longitudeName='lon',
 latitudeName='lat',
 hasVerticalDimension=.true.,
 verticalName='lev',
/

&forecast
 numberOfTimes=1,
 validTimes(1)='${timeStamp1}',
 timeIndices(1) = 1,
 forecastHours(1)=${forecastHour},

/

&variables
 numberOfVariables=8,

 variableRanks(1) = 4,
 variableLevelTypes(1) = 11,
 variableNamesIn(1)='PL',
 variableNamesOut(1)='PRESSURE',
 variableUnits(1)='Pa',
 variableDescriptions(1)='Pressure',

 variableRanks(2) = 4,
 variableLevelTypes(2) = 11,
 variableNamesIn(2)='T',
 variableNamesOut(2)='TT',
 variableUnits(2)='K',
 variableDescriptions(2)='Temperature',

 variableRanks(3) = 4,
 variableLevelTypes(3) = 11,
 variableNamesIn(3)='U',
 variableNamesOut(3)='UU',
 variableUnits(3)='m s-1',
 variableDescriptions(3)='U',

 variableRanks(4) = 4,
 variableLevelTypes(4) = 11,
 variableNamesIn(4)='V',
 variableNamesOut(4)='VV',
 variableUnits(4)='m s-1',
 variableDescriptions(4)='V',

 variableRanks(5) = 4,
 variableLevelTypes(5) = 11,
 variableNamesIn(5)='QV',
 variableNamesOut(5)='SPECHUMD',
 variableUnits(5)='kg kg-1',
 variableDescriptions(5)='Specific Humidity',

 variableRanks(6) = 4,
 variableLevelTypes(6) = 11,
 variableNamesIn(6)='H',
 variableNamesOut(6)='HGT',
 variableUnits(6)='m',
 variableDescriptions(6)='Height',

 variableRanks(7) = 3,
 variableLevelTypes(7) = 4,
 variableNamesIn(7)='SLP',
 variableNamesOut(7)='PMSL',
 variableUnits(7)='Pa',
 variableDescriptions(7)='Sea-level Pressure',

 variableRanks(8) = 3,
 variableLevelTypes(8) = 1,
 variableNamesIn(8)='PS',
 variableNamesOut(8)='PSFC',
 variableUnits(8)='Pa',
 variableDescriptions(8)='Surface Pressure',

/

EOF

    $converter
    if [ $? -ne 0 ] ; then
	echo "ERROR returned by $converter"
	exit 1
    fi

    # Go to next time
    currentTimeString=`${updateDateTime} $currentTimeString $deltaHour`
done

exit 0

