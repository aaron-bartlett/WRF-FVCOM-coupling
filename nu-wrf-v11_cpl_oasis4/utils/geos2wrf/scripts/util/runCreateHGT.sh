#!/bin/sh
#
# Sample script for running createHGT

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

    cat DELP_MODEL_LEVEL:${timeStamp} \
        TT_MODEL_LEVEL:${timeStamp} \
        SPECHUMD_MODEL_LEVEL:${timeStamp} \
        SOILHGT_GROUND_LEVEL:${timeStamp} \
        > GEOS_TMP:${timeStamp}

    cat > namelist.createHGT<<EOF
&input
  directory='./',
  prefix='GEOS_TMP',
  year=${cyear}, 
  month=${cmonth},
  day=${cday},
  hour=${chour},
  layerPressureThicknessName='DELP',
  layerTemperatureName='TT',
  layerSpecificHumidityName='SPECHUMD',
  soilHeightName='SOILHGT',
  modelTopPressure=1.,
/

EOF

    $height
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from $height"
	exit 1
    fi
#    mv HGT:${timeStamp} HGT_MODEL_LEVEL:${timeStamp}
    rm GEOS_TMP:${timeStamp}

    # Go to next time
    currentTimeString=`${updateDateTime} $currentTimeString $deltaHour`
done

exit 0

