#!/bin/sh
#
# Sample script for running createLANDSEA

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

    cat FRLAKE_GROUND_LEVEL:${timeStamp} \
        FROCEAN_GROUND_LEVEL:${timeStamp} \
        > GEOS_TMP:${timeStamp}

    cat > namelist.createLANDSEA<<EOF
&input
  directory='./',
  prefix='GEOS_TMP',
  year=${cyear}, 
  month=${cmonth},
  day=${cday},
  hour=${chour},
  lakeFractionName='FRLAKE',
  oceanFractionName='FROCEAN',
/

EOF

    $landsea
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from $landsea"
	exit 1
    fi
#    mv LANDSEA:${timeStamp} LANDSEA_GROUND_LEVEL:${timeStamp}
    rm GEOS_TMP:${timeStamp}

    # Go to next time
    currentTimeString=`${updateDateTime} $currentTimeString $deltaHour`
done

exit 0

