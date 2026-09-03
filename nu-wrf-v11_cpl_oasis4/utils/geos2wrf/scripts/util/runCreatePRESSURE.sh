#!/bin/sh
#
# Sample script for running createPRESSURE

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

    cp DELP_MODEL_LEVEL:${timeStamp} GEOS_TMP:${timeStamp}

    cat > namelist.createPRESSURE<<EOF
&input
  directory='./',
  prefix='GEOS_TMP',
  year=${cyear}, 
  month=${cmonth},
  day=${cday},
  hour=${chour},
  layerPressureThicknessName='DELP',
  modelTopPressure=1.,
/

EOF

    $pressure
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from $pressure"
	exit 1
    fi
#    mv PRESSURE:${timeStamp} PRESSURE_MODEL_LEVEL:${timeStamp}
    rm GEOS_TMP:${timeStamp}

    # Go to next time
    currentTimeString=`${updateDateTime} $currentTimeString $deltaHour`
done

exit 0

