#!/bin/sh 
# Sample script for processing MERRA inst6_3d_ana_Np file with geos2wps,
# Assumes MERRA file has data for 00Z, 06Z, 12Z, and 18Z.
# NOTE:  This version processes the isobaric (constant pressure) levels.

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
  geosFileFormat=2,
  geosFileName='MERRA300.prod.assim.inst6_3d_ana_Np.${time}.SUB.nc.1'
  outputDirectory='./',
/

&coordinates
 longitudeName='longitude',
 latitudeName='latitude',
 hasVerticalDimension=.true.,
 verticalName='levels',
/

&forecast
 numberOfTimes=1,
 validTimes(1)='${timeStamp1}',
 timeIndices(1) = ${timeLevel},,
 forecastHours(1)=0,
/

&variables
 numberOfVariables=7,

 variableRanks(1) = 3,
 variableLevelTypes(1) = 4,
 variableNamesIn(1)='slp',
 variableNamesOut(1)='PMSL',
 variableUnits(1)='Pa',
 variableDescriptions(1)='Sea-level Pressure',

 variableRanks(2) = 3,
 variableLevelTypes(2) = 1,
 variableNamesIn(2)='ps',
 variableNamesOut(2)='PSFC',
 variableUnits(2)='Pa',
 variableDescriptions(2)='Surface Pressure',

 variableRanks(3) = 4,
 variableLevelTypes(3) = 12,
 variableNamesIn(3)='h',
 variableNamesOut(3)='HGT',
 variableUnits(3)='m',
 variableDescriptions(3)='Height',

 variableRanks(4) = 4,
 variableLevelTypes(4) = 12,
 variableNamesIn(4)='t',
 variableNamesOut(4)='TT',
 variableUnits(4)='K',
 variableDescriptions(4)='Temperature',

 variableRanks(5) = 4,
 variableLevelTypes(5) = 12,
 variableNamesIn(5)='qv',
 variableNamesOut(5)='SPECHUMD',
 variableUnits(5)='kg kg-1',
 variableDescriptions(5)='Specific Humidity',

 variableRanks(6) = 4,
 variableLevelTypes(6) = 12,
 variableNamesIn(6)='u',
 variableNamesOut(6)='UU',
 variableUnits(6)='m s-1',
 variableDescriptions(6)='U',

 variableRanks(7) = 4,
 variableLevelTypes(7) = 12,
 variableNamesIn(7)='v',
 variableNamesOut(7)='VV',
 variableUnits(7)='m s-1',
 variableDescriptions(7)='V',

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

