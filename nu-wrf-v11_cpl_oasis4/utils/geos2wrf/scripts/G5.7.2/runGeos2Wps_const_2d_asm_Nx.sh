#!/bin/sh 
#
# Sample script for processing GEOS const_2d_asm_Nx file with geos2wps,

startTimeString=$1
endTimeString=$2
deltaHour=$3
cfgFile=$4

source $cfgFile

currentTimeString=${startTimeString}
while [ "$currentTimeString" -le "$endTimeString" ] ; do

    cyear=`echo $currentTimeString | cut -c1-4`
    cmonth=`echo $currentTimeString | cut -c5-6`
    cday=`echo $currentTimeString | cut -c7-8`
    chour=`echo $currentTimeString | cut -c9-10`
    timeStamp="${cyear}-${cmonth}-${cday}_${chour}"

    currentTimeString="${cyear}${cmonth}${cday}${chour}"
    forecastHour=`${calcForecastHour} $startTimeString $currentTimeString`
 
    cat > namelist.geos2wps<<EOF
&files
  geosFileFormat=2,
  geosFileName='DAS.ops.asm.const_2d_asm_Nx.GEOS572.00000000_0000.V01.nc4'
  outputDirectory='./',
/

&coordinates
 longitudeName='lon',
 latitudeName='lat',
/

&forecast
 numberOfTimes=1,
 validTimes(1)='${timeStamp}',
 timeIndices(1) = 1,
 forecastHours(1)=${forecastHour},
/

&variables
 numberOfVariables=3,

 variableRanks(1) = 3,
 variableLevelTypes(1) = 1,
 variableNamesIn(1)='PHIS',
 variableNamesOut(1)='PHIS',
 variableUnits(1)='m**2 s**-2',
 variableDescriptions(1)='Surface geopotential',

 variableRanks(2) = 3,
 variableLevelTypes(2) = 1,
 variableNamesIn(2)='FRLAKE',
 variableNamesOut(2)='FRLAKE',
 variableUnits(2)='proprtn',
 variableDescriptions(2)='Lake fraction',

 variableRanks(3) = 3,
 variableLevelTypes(3) = 1,
 variableNamesIn(3)='FROCEAN',
 variableNamesOut(3)='FROCEAN',
 variableUnits(3)='proprtn',
 variableDescriptions(3)='Ocean fraction',

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

