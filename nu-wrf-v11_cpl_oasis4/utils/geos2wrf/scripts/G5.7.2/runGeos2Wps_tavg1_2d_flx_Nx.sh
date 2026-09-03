#!/bin/sh 
# Sample script for processing GEOS tavg1_2d_flx_Nx file with geos2wps.

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
lastOne=0
while [ "$currentTimeString" -le "$endTimeString" ] ; do

    if [ "$currentTimeString" -eq "$endTimeString" ] ; then

	nominalCYear=`echo $currentTimeString | cut -c1-4`
	nominalCMonth=`echo $currentTimeString | cut -c5-6`
	nominalCDay=`echo $currentTimeString | cut -c7-8`
	nominalCHour=`echo $currentTimeString | cut -c9-10`

	fullForecastLength=`$calcForecastHour $startTimeString $endTimeString`
	forecastLength=$((fullForecastLength-1))
	currentTimeString=`${updateDateTime} $startTimeString $forecastLength`
	lastOne=1
    fi

    cyear=`echo $currentTimeString | cut -c1-4`
    cmonth=`echo $currentTimeString | cut -c5-6`
    cday=`echo $currentTimeString | cut -c7-8`
    chour=`echo $currentTimeString | cut -c9-10`
    
    if [ $lastOne -eq 0 ] ; then
	nominalCYear=$cyear
	nominalCMonth=$cmonth
	nominalCDay=$cday
	nominalCHour=$chour
    fi

    timeStamp="${nominalCYear}-${nominalCMonth}-${nominalCDay}_${nominalCHour}"

    cdate="${cyear}${cmonth}${cday}"

    timeStamp1=$timeStamp

    currentTimeString="${cyear}${cmonth}${cday}${chour}"
    if [ $lastOne -eq 0 ] ; then
	forecastHour=`${calcForecastHour} $startTimeString $currentTimeString`
    else
	forecastHour=${fullForecastLength}
    fi

    cat > namelist.geos2wps<<EOF
&files
  geosFileFormat=2,
  geosFileName='e572p5_fp.tavg1_2d_flx_Nx.${sdate}_${shour}z+${cdate}_${chour}30z.nc4'
  outputDirectory='./',
/

&coordinates
 longitudeName='lon',
 latitudeName='lat',
/

&forecast
 numberOfTimes=1,
 validTimes(1)='${timeStamp1}',
 timeIndices(1) = 1,
 forecastHours(1)=${forecastHour},
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

    if [ $lastOne -eq 1 ] ; then
	break
    fi
	
    # Go to next time
    currentTimeString=`${updateDateTime} $currentTimeString $deltaHour`
done

exit 0

