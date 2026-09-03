#!/bin/sh -x

# Sample script for downloading operational GEOS G5.7.2 data via ftp

#startDateTime="2012051400"
startDateTime=$1
cfgFile=$2

source $cfgFile

sYear=`echo $startDateTime | cut -c1-4`
sMonth=`echo $startDateTime | cut -c5-6`
sDay=`echo $startDateTime | cut -c7-8`
sHour=`echo $startDateTime | cut -c9-10`

sDate=${sYear}${sMonth}${sDay}

# First, get const_2d_asm_Nx file
top=ftp://gmao_ops@ftp.nccs.nasa.gov/fp/das/
wget ${top}/DAS.ops.asm.const_2d_asm_Nx.GEOS572.00000000_0000.V01.nc4 \
     --password=""
if [ $? -ne 0 ] ; then
    echo "ERROR returned from wget!"
    exit 1
fi


# Next, get all model level data files (inst3_3d_asm_Nv) for current run
top=ftp://gmao_ops@ftp.nccs.nasa.gov/fp/forecast/
prefix=e572p5_fp
col=inst3_3d_asm_Nv
wget "$top/Y${sYear}/M${sMonth}/D${sDay}/H${sHour}/${prefix}.${col}.${sDate}_${sHour}z*.nc4" --password=""
if [ $? -ne 0 ] ; then
    echo "ERROR returned from wget!"
    exit 1
fi

# Next, get tavg1 files.  We only want those files 30 minutes after
# each inst3_3d_asm_Nv valid time, except we also want the very last dump.

top=ftp://gmao_ops@ftp.nccs.nasa.gov/fp/forecast/
prefix=e572p5_fp

for col in tavg1_2d_flx_Nx tavg1_2d_slv_Nx ; do

    for validTime in 0030 0330 0630 0930 1230 1530 1830 2130 ; do
	wget "$top/Y${sYear}/M${sMonth}/D${sDay}/H${sHour}/${prefix}.${col}.${sDate}_${sHour}z*${validTime}z.nc4" --password=""
	if [ $? -ne 0 ] ; then
	    echo "ERROR returned from wget!"
	    exit 1
	fi
    done

    if [ "$sHour" -eq 0 ] || [ "$sHour" -eq 12 ] ; then
	lastDate=`${updateDateTime} $startDateTime 119`
    else
	lastDate=`${updateDateTime} $startDateTime 29`
    fi

    lYear=`echo $lastDate | cut -c1-4`
    lMonth=`echo $lastDate | cut -c5-6`
    lDay=`echo $lastDate | cut -c7-8`
    lHour=`echo $lastDate | cut -c9-10`

    lDate=${lYear}${lMonth}${lDay}

    wget $top/Y${sYear}/M${sMonth}/D${sDay}/H${sHour}/${prefix}.${col}.\
${sDate}_${sHour}z+${lDate}_${lHour}30z.nc4 --password=""
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from wget!"
	exit 1
    fi

done

exit 0