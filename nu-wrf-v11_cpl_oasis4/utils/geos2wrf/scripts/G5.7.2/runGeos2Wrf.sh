#!/bin/sh
# Sample driver script for processing GEOS files with geos2wrf programs.

#startTimeString="2012051400"
#endTimeString="2012051900"
#deltaHour="03"

#source ./geos2wrf.cfg

startTimeString=$1
endTimeString=$2
deltaHour=$3
cfgFile=$4

source $cfgFile

${const_2d_asm_Nx} $startTimeString $endTimeString $deltaHour $cfgFile
if [ $? -ne 0 ] ; then
   echo "ERROR processing const_2d_asm_Nx files!"
fi

${inst6_3d_ana_Nv} $startTimeString $endTimeString $deltaHour $cfgFile
if [ $? -ne 0 ] ; then
    echo "ERROR processing inst6_3d_ana_Nv files!"
fi

${tavg1_2d_flx_Nx} $startTimeString $endTimeString $deltaHour $cfgFile
if [ $? -ne 0 ] ; then
    echo "ERROR processing tavg1_2d_slv_Nx files!"
fi

${tavg1_2d_slv_Nx} $startTimeString $endTimeString $deltaHour $cfgFile
if [ $? -ne 0 ] ; then
    echo "ERROR processing tavg1_2d_slv_Nx files!"
fi

${runCreateSOILHGT} $startTimeString $endTimeString $deltaHour $cfgFile
if [ $? -ne 0 ] ; then
    echo "ERROR generating SOILHGT files!"
fi

${runCreateLANDSEA} $startTimeString $endTimeString $deltaHour $cfgFile
if [ $? -ne 0 ] ; then
   echo "ERROR generating LANDSEA files!"
fi

${runCreateRH} $startTimeString $endTimeString $deltaHour $cfgFile
if [ $? -ne 0 ] ; then
   echo "ERROR generating RH files!"
fi

currentTimeString=${startTimeString}
while [ "$currentTimeString" -le "$endTimeString" ] ; do

    cyear=`echo $currentTimeString | cut -c1-4`
    cmonth=`echo $currentTimeString | cut -c5-6`
    cday=`echo $currentTimeString | cut -c7-8`
    chour=`echo $currentTimeString | cut -c9-10`
    timeStamp="${cyear}-${cmonth}-${cday}_${chour}"

    cat HGT_MODEL_LEVEL:${timeStamp} \
        LANDSEA_GROUND_LEVEL:${timeStamp} \
	PMSL_MEAN_SEA_LEVEL:${timeStamp} \
	PRESSURE_MODEL_LEVEL:${timeStamp} \
	PSFC_GROUND_LEVEL:${timeStamp} \
	RH_MODEL_LEVEL:${timeStamp} \
	RH_2M_ABOVE_GROUND_LEVEL:${timeStamp} \
	SEAICE_GROUND_LEVEL:${timeStamp} \
	SKINTEMP_GROUND_LEVEL:${timeStamp} \
	SOILHGT_GROUND_LEVEL:${timeStamp} \
	TT_2M_ABOVE_GROUND_LEVEL:${timeStamp} \
	TT_MODEL_LEVEL:${timeStamp} \
	UU_10M_ABOVE_GROUND_LEVEL:${timeStamp} \
	UU_MODEL_LEVEL:${timeStamp} \
	VV_10M_ABOVE_GROUND_LEVEL:${timeStamp} \
	VV_MODEL_LEVEL:${timeStamp} \
	> GEOS:${timeStamp}.save

    rm *:${timeStamp}
    mv GEOS:${timeStamp}.save GEOS:${timeStamp}

    # Go to next time
    currentTimeString=`${updateDateTime} $currentTimeString $deltaHour`
done

exit 0
