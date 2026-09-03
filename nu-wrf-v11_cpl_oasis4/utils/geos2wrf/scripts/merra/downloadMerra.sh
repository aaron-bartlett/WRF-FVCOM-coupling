#!/bin/sh

# Sample script for downloading MERRA data.  Very crude; user will likely need
# to tweak server name and file prefix depending on collection and valid dates.

startDate=$1
endDate=$2
cfgFile=$3

source $cfgFile

# First, get const_2d_asm_Nx file
top=ftp://goldsmr2.sci.gsfc.nasa.gov/data/s4pa/MERRA_MONTHLY/MAC0NXASM.5.2.0
wget ${top}/1979/MERRA300.prod.assim.const_2d_asm_Nx.00000000.hdf
if [ $? -ne 0 ] ; then
    echo "ERROR returned from wget!"
    exit 1
fi

# Next, get all model level data files (inst6_3d_asm_Nv)
currentDate=$startDate
top=http://goldsmr3.sci.gsfc.nasa.gov/opendap/MERRA/MAI6NVANA.5.2.0
while [ "$currentDate" -le "$endDate" ] ; do
    
    cyear=`echo $currentDate | cut -c1-4`
    cmonth=`echo $currentDate | cut -c5-6`
    cday=`echo $currentDate | cut -c7-8`

    wget ${top}/${cyear}/${cmonth}/MERRA300.prod.assim.inst6_3d_ana_Nv.${cyear}${cmonth}${cday}.hdf
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from wget!"
	exit 1
    fi

    # Go to next day
    currentDate=`${updateDateTime} ${currentDate}00 24`
    currentDate=`echo $currentDate | cut -c1-8`
done

# Next, get all model level data files (inst6_3d_asm_Np)
currentDate=$startDate
top=http://goldsmr3.sci.gsfc.nasa.gov/opendap/MERRA/MAI6NPANA.5.2.0
while [ "$currentDate" -le "$endDate" ] ; do
    
    cyear=`echo $currentDate | cut -c1-4`
    cmonth=`echo $currentDate | cut -c5-6`
    cday=`echo $currentDate | cut -c7-8`

    wget ${top}/${cyear}/${cmonth}/MERRA300.prod.assim.inst6_3d_ana_Np.${cyear}${cmonth}${cday}.hdf
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from wget!"
	exit 1
    fi

    # Go to next day
    currentDate=`${updateDateTime} ${currentDate}00 24`
    currentDate=`echo $currentDate | cut -c1-8`
done

# Next, get all model level data files (tavg1_2d_ocn_Nx)
currentDate=$startDate
top=http://goldsmr2.sci.gsfc.nasa.gov/opendap/MERRA/MAT1NXOCN.5.2.0
while [ "$currentDate" -le "$endDate" ] ; do
    
    cyear=`echo $currentDate | cut -c1-4`
    cmonth=`echo $currentDate | cut -c5-6`
    cday=`echo $currentDate | cut -c7-8`

    wget ${top}/${cyear}/${cmonth}/MERRA300.prod.assim.tavg1_2d_ocn_Nx.${cyear}${cmonth}${cday}.hdf
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from wget!"
	exit 1
    fi

    # Go to next day
    currentDate=`${updateDateTime} ${currentDate}00 24`
    currentDate=`echo $currentDate | cut -c1-8`
done

# Next, get all model level data files (tavg1_2d_slv_Nx)
currentDate=$startDate
top=http://goldsmr2.sci.gsfc.nasa.gov/opendap/MERRA/MAT1NXSLV.5.2.0
while [ "$currentDate" -le "$endDate" ] ; do
    
    cyear=`echo $currentDate | cut -c1-4`
    cmonth=`echo $currentDate | cut -c5-6`
    cday=`echo $currentDate | cut -c7-8`

    wget ${top}/${cyear}/${cmonth}/MERRA300.prod.assim.tavg1_2d_slv_Nx.${cyear}${cmonth}${cday}.hdf
    if [ $? -ne 0 ] ; then
	echo "ERROR returned from wget!"
	exit 1
    fi

    # Go to next day
    currentDate=`${updateDateTime} ${currentDate}00 24`
    currentDate=`echo $currentDate | cut -c1-8`
done
