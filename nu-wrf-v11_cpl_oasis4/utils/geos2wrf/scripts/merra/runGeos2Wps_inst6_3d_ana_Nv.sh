#!/bin/sh 
# Sample script for processing MERRA inst6_3d_ana_Nv file with geos2wps.
# Assumes MERRA file has data for 00Z, 06Z, 12Z, and 18Z.

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

    if [ "$chour" -ne 0 ] && [ "$chour" -ne 6 ] && [ "$chour" -ne 12 ] && [ "$chour" -ne 18 ] ; then
	echo "ERROR, invalid hour for MERRA data!"
	echo "hour = $chour"
	echo "Must be 00, 06, 12, or 18!"
	exit 1
    fi

    timeStamp="${cyear}-${cmonth}-${cday}_${chour}"

    time="${cyear}${cmonth}${cday}"

    timeStamp1=$timeStamp
    timeLevel=$((chour/6+1))

    cat > namelist.geos2wps<<EOF
&files
  geosFileFormat=4,
  geosFileName='MERRA300.prod.assim.inst6_3d_ana_Nv.${time}.hdf'
  outputDirectory='./',
/

&coordinates
 longitudeName='XDim',
 latitudeName='YDim',
 hasVerticalDimension=.true.,
 verticalName='Height',
/

&forecast
 numberOfTimes=1,
 validTimes(1)='${timeStamp1}',
 timeIndices(1) = ${timeLevel},
 forecastHours(1)=0,
/

&variables
 numberOfVariables=6,

 variableRanks(1) = 4,
 variableLevelTypes(1) = 11,
 variableNamesIn(1)='DELP',
 variableNamesOut(1)='DELP',
 variableUnits(1)='Pa',
 variableDescriptions(1)='Layer pressure thickness',

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

 variableRanks(6) = 3,
 variableLevelTypes(6) = 1,
 variableNamesIn(6)='PS',
 variableNamesOut(6)='PSFC',
 variableUnits(6)='Pa',
 variableDescriptions(6)='Surface Pressure',

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

