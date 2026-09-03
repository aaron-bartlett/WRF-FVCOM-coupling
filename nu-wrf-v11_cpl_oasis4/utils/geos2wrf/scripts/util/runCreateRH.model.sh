#!/bin/sh
#
# Sample script for running createRH
# NOTE:  This processes upper air data on GEOS model levels.

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

    cat PSFC_GROUND_LEVEL:${timeStamp} \
        PRESSURE_MODEL_LEVEL:${timeStamp} \
        TT_MODEL_LEVEL:${timeStamp} \
	TT_2M_ABOVE_GROUND_LEVEL:${timeStamp} \
        SPECHUMD_MODEL_LEVEL:${timeStamp} \
	SPECHUMD_2M_ABOVE_GROUND_LEVEL:${timeStamp} \
        > GEOS_TMP:${timeStamp}

    cat > namelist.createRH<<EOF
&input
  directory='./',
  prefix='GEOS_TMP',
  year=${cyear}, 
  month=${cmonth},
  day=${cday},
  hour=${chour},
  processSurfacePressure=.true.,
  surfacePressureName='PSFC',
  pressureName='PRESSURE',
  temperatureName='TT',
  specificHumidityName='SPECHUMD',
/

EOF

    $rh
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from $rh"
	exit 1
    fi
    rm GEOS_TMP:${timeStamp}

    # Go to next time
    currentTimeString=`${updateDateTime} $currentTimeString $deltaHour`
done

exit 0

