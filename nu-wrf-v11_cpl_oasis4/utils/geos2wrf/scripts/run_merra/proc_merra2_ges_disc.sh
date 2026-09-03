#!/bin/sh
#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Science and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# SCRIPT: proc_merra2_ges_disc.sh
#
# AUTHOR:
# Eric Kemp, NASA CISTO/SSAI
#
# DESCRIPTION:
# Script for fetching MERRA2 data from NASA GES DISC web site for use by
# MERRA2WRF. Requires wget utility and Get_dates_daily.py script. Based on 
# Jossy Jacob's Run_MERRA2.csh script.
#
# REVISION HISTORY:
# 18 Sep 2015 - First version.
#
#------------------------------------------------------------------------------

# Process command line arguments.
if [ "$#" -ne 4 ] ; then
    echo 'Usage: proc_merra2_ges_disc.sh STARTDATE ENDDATE RUNDIR NUWRFDIR'
    echo 'Example: ./proc_merra2_ges_disc.sh 19890101 19890104 MERRA_data /path/to/nuwrf'
    exit 1
fi
STARTDATE=$1
ENDDATE=$2
RUNDIR=$3
NUWRFDIR=$4

# Set WORKDIR and ensure it is an absolute path. Create directory and 
# subdirectories as needed.
WORKDIR="$RUNDIR"
if [ ! -e $WORKDIR ] ; then
    mkdir $WORKDIR || exit 1
fi
cd $WORKDIR
WORKDIR=`pwd` # Absolute path.
if [ ! -e $WORKDIR/m2wIn ] ; then
    mkdir $WORKDIR/m2wIn || exit 1
fi
if [ ! -e $WORKDIR/m2wOut ] ; then
    mkdir $WORKDIR/m2wOut || exit 1
fi

# Link executables from NUWRFDIR
ln -fs $NUWRFDIR/utils/bin/merra2wrf.x $WORKDIR/merra2wrf.x || exit 1
if [ ! -e $WORKDIR/merra2wrf.x ] ; then
    echo "ERROR, merra2wrf.x does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/utils/merra2wrf/scripts/run_merra/Get_dates_daily.py \
   $WORKDIR/Get_dates_daily.py || exit 1
if [ ! -e $WORKDIR/Get_dates_daily.py ] ; then
    echo "ERROR, $WORKDIR/Get_dates_daily.py does not exist!"
    exit 1
fi

# Construct list of dates.
cd $WORKDIR || exit 1
DT=1 # Day
$WORKDIR/Get_dates_daily.py $STARTDATE $ENDDATE $DT || exit 1
if [ ! -e "$WORKDIR/Datetime.dat" ] ; then
    echo "ERROR, Datetime.dat not generated!"
    exit 1
fi

# Get const_2d_asm_Nx file. This is used for all dates.
cd $WORKDIR/m2wIn || exit 1
# Updated 21 March 2016
#file1=MERRA2_100.const_2d_asm_Nx.00000000.nc4 
file1=MERRA2_101.const_2d_asm_Nx.00000000.nc4 
if [ ! -e "$file1" ] ; then
    WEBSITE=ftp://goldsmr4.sci.gsfc.nasa.gov
    path=data/s4pa/MERRA2_MONTHLY/M2C0NXASM.5.12.4/1980/$file1
    wget $WEBSITE/$path || exit 1
fi

# Loop through the dates, fetch date specific files, and run MERRA2WRF.
for date in $(cat $WORKDIR/Datetime.dat) ; do

    YEAR=`echo $date | cut -c1-4`
    MONTH=`echo $date | cut -c5-6`
    DAY=`echo $date | cut -c7-8`

    # File prefix depends on MERRA2 "data stream" which depends on year.
    if [ "$YEAR" -lt 1992 ] ; then
	MERRANAME="MERRA2_100"
    elif [ "$YEAR" -lt 2001 ] ; then
	MERRANAME="MERRA2_200"
    elif [ "$YEAR" -lt 2011 ] ; then
	MERRANAME="MERRA2_300"
    else
	MERRANAME="MERRA2_400"
    fi

    # Get inst6_3d_ana_Nv file.
    file2=$MERRANAME.inst6_3d_ana_Nv.$YEAR$MONTH$DAY.nc4
    if [ ! -e "$file2" ] ; then
	WEBSITE=ftp://goldsmr5.sci.gsfc.nasa.gov
	path=data/s4pa/MERRA2/M2I6NVANA.5.12.4/$YEAR/$MONTH/$file2
	wget $WEBSITE/$path || exit 1
    fi

    # Get inst6_3d_ana_Np file.
    file3=$MERRANAME.inst6_3d_ana_Np.$YEAR$MONTH$DAY.nc4
    if [ ! -e "$file3" ] ; then
	WEBSITE=ftp://goldsmr5.sci.gsfc.nasa.gov
	path=data/s4pa/MERRA2/M2I6NPANA.5.12.4/$YEAR/$MONTH/$file3
	wget $WEBSITE/$path || exit 1
    fi

    # Get tavg1_2d_slv_Nx file.
    file4=$MERRANAME.tavg1_2d_slv_Nx.$YEAR$MONTH$DAY.nc4
    if [ ! -e "$file4" ] ; then
	WEBSITE=ftp://goldsmr4.sci.gsfc.nasa.gov
	path=data/s4pa/MERRA2/M2T1NXSLV.5.12.4/$YEAR/$MONTH/$file4
	wget $WEBSITE/$path || exit 1
    fi

    # Get tavg1_2d_ocn_Nx file.
    file5=$MERRANAME.tavg1_2d_ocn_Nx.$YEAR$MONTH$DAY.nc4
    if [ ! -e "$file5" ] ; then
	WEBSITE=ftp://goldsmr4.sci.gsfc.nasa.gov
	path=data/s4pa/MERRA2/M2T1NXOCN.5.12.4/$YEAR/$MONTH/$file5
	wget $WEBSITE/$path || exit 1
    fi
    
    MERRADATE="$YEAR-$MONTH-$DAY"
    MERRADATE2="$YEAR$MONTH$DAY"

    # Now run MERRA2WRF for the collected data.
    cd $WORKDIR || exit 1
    cat > namelist.merra2wrf_daily <<EOF
&input
    outputDirectory = '$WORKDIR/m2wOut',
    merraDirectory = '$WORKDIR/m2wIn',
    merraFormat_const_2d_asm_Nx = 2,
    merraFile_const_2d_asm_Nx = '$file1',
    numberOfDays=1,
    merraDates(1)="$MERRADATE",
    merraFormat_inst6_3d_ana_Nv = 2,
    merraFiles_inst6_3d_ana_Nv(1) = '$file2',
    merraFormat_inst6_3d_ana_Np = 2,
    merraFiles_inst6_3d_ana_Np(1) = '$file3',
    merraFormat_tavg1_2d_slv_Nx = 2,
    merraFiles_tavg1_2d_slv_Nx(1) = '$file4',
    merraFormat_tavg1_2d_ocn_Nx = 2,
    merraFiles_tavg1_2d_ocn_Nx(1) = '$file5',
/
EOF
    $WORKDIR/merra2wrf namelist.merra2wrf_daily || exit 1

    # Daily file cleanup.
    rm $WORKDIR/m2wIn/$file2
    rm $WORKDIR/m2wIn/$file3
    rm $WORKDIR/m2wIn/$file4
    rm $WORKDIR/m2wIn/$file5

    cd $WORKDIR/m2wIn || exit 1

done

# Final file cleanup.
rm $WORKDIR/m2wIn/$file1

# The End.
echo "MERRA2WRF output files are in $WORKDIR/m2wOut"
echo "Completed MERRA2 processing"
exit 0


