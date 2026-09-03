#!/bin/sh 
# Sample script for processing MERRA tavg1_2d_slv_Nx file with geos2wps,
# Assumes MERRA file has hourly data from 00:30 to 23:30Z.

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
    timeLevel=$((chour+1)) # Hourly data in this file.

    cat > namelist.geos2wps<<EOF
&files
  geosFileFormat=4,
  geosFileName='MERRA300.prod.assim.tavg1_2d_slv_Nx.${time}.hdf'
  outputDirectory='./',
/

&coordinates
 longitudeName='XDim',
 latitudeName='YDim',
/

&forecast
 numberOfTimes=1,
 validTimes(1)='${timeStamp1}',
 timeIndices(1) = ${timeLevel},
 forecastHours(1)=0,
/

&variables
 numberOfVariables=5,

 variableRanks(1) = 3,
 variableLevelTypes(1) = 3,
 variableNamesIn(1)='U10M',
 variableNamesOut(1)='UU',
 variableUnits(1)='m s-1',
 variableDescriptions(1)='U                 at 10 m',

 variableRanks(2) = 3,
 variableLevelTypes(2) = 3,
 variableNamesIn(2)='V10M',
 variableNamesOut(2)='VV',
 variableUnits(2)='m s-1',
 variableDescriptions(2)='V                 at 10 m',

 variableRanks(3) = 3,
 variableLevelTypes(3) = 2,
 variableNamesIn(3)='T2M',
 variableNamesOut(3)='TT',
 variableUnits(3)='K',
 variableDescriptions(3)='Temperature       at 2 m',

 variableRanks(4) = 3,
 variableLevelTypes(4) = 2,
 variableNamesIn(4)='QV2M',
 variableNamesOut(4)='SPECHUMD',
 variableUnits(4)='kg kg-1',
 variableDescriptions(4)='Specific Humidity at 2 m',

 variableRanks(5) = 3,
 variableLevelTypes(5) = 1,
 variableNamesIn(5)='TS',
 variableNamesOut(5)='SKINTEMP',
 variableUnits(5)='K',
 variableDescriptions(5)='Skin temperature (can use for SST also)',
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

