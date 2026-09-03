#!/bin/sh 
# Sample script for processing MERRA tavg1_2d_ocn_Nx file with geos2wps,
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
  geosFileName='MERRA300.prod.assim.tavg1_2d_ocn_Nx.${time}.hdf'
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
 numberOfVariables=1,

 variableRanks(1) = 3,
 variableLevelTypes(1) = 1,
 variableNamesIn(1)='FRSEAICE',
 variableNamesOut(1)='SEAICE',
 variableUnits(1)='proprtn',
 variableDescriptions(1)='Ice flag',
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

