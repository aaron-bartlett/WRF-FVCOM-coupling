#!/bin/sh
#
# Sample script for running extrapIsobaric

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

    cat HGT_ISOBARIC_LEVEL:${timeStamp} \
        TT_ISOBARIC_LEVEL:${timeStamp} \
        RH_ISOBARIC_LEVEL:${timeStamp} \
        UU_ISOBARIC_LEVEL:${timeStamp} \
        VV_ISOBARIC_LEVEL:${timeStamp} \
        > GEOS_TMP:${timeStamp}

    cat > namelist.extrapIsobaric<<EOF
&input
  directory='./',
  prefix='GEOS_TMP',
  year=${cyear}, 
  month=${cmonth},
  day=${cday},
  hour=${chour},
  geopotentialHeightName='HGT',
  temperatureName='TT',
  relativeHumidityName='RH',
  uName='UU',
  vName='VV',
/

EOF

    $extrap
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from $rh"
	exit 1
    fi
    rm GEOS_TMP:${timeStamp}

    # Go to next time
    currentTimeString=`${updateDateTime} $currentTimeString $deltaHour`
done

exit 0

