#!/bin/sh
#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Science and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#
# SCRIPT: proc_merra2_gocart_ges_disc.sh
#
# AUTHOR:
# Eric Kemp, NASA CISTO/SSAI
#
# DESCRIPTION:
# Script for fetching MERRA2 GOCART data from NASA GES DISC web site and
# processing with GOCART2WRF. Requires wget utility and Get_dates_daily.py 
# script. Based on Jossy Jacob's Run_MERRA2.csh script.
#
# REVISION HISTORY:
# 21 Sep 2015 - First version.
#
#------------------------------------------------------------------------------

# Process command line arguments.
if [ "$#" -ne 5 ] ; then
    echo 'Usage: proc_merra2_gocart_ges_disc.sh NUMDOMAINS STARTDATE ENDDATE RUNDIR NUWRFDIR'
    echo 'Example: ./proc_merra2_gocart_ges_disc.sh 2 19890101 19890104 MERRA_data /path/to/nuwrf'
    exit 1
fi
MAXDOM=$1
STARTDATE=$2
ENDDATE=$3
RUNDIR=$4
NUWRFDIR=$5

# Set WORKDIR and ensure it is an absolute path. Create directory and 
# subdirectories as needed.
WORKDIR="$RUNDIR"
if [ ! -e $WORKDIR ] ; then
    mkdir $WORKDIR || exit 1
fi
cd $WORKDIR
WORKDIR=`pwd` # Absolute path.

# Make sure wrfinput and wrfbdy files are in present directory
cd $WORKDIR
if [ ! -e ./wrfbdy_d01 ] ; then
    echo "ERROR, wrfbdy_d01 is not in $WORKDIR!"
    exit 1
fi
cp ./wrfbdy_d01 wrfbdy_d01.real || exit 1
for (( x=1; x <= $MAXDOM; x += 1 )) ; do
    printf -v domain %2.2d $x
    if [ ! -e ./wrfinput_d$domain ] ; then
	echo "ERROR, wrinput_d$domain is not in $WORKDIR!"
	exit 1
    fi
    cp ./wrfinput_d$domain ./wrfinput_d$domain.real || exit 1
done

# Link executables from NUWRFDIR
if [ ! -e $NUWRFDIR/utils/bin/gocart2wrf.x ] ; then
    echo $NUWRFDIR/utils/bin/gocart2wrf.x does not exist
    exit 1   
fi
ln -fs $NUWRFDIR/utils/bin/gocart2wrf.x || exit 1

if [ ! -e $NUWRFDIR/utils/geos2wrf/scripts/run_merra/Get_dates_daily.py ] ; then
    echo $NUWRFDIR/utils/geos2wrf/scripts/run_merra/Get_dates_daily.py does not exist
    exit 1   
fi
ln -fs $NUWRFDIR/utils/geos2wrf/scripts/run_merra/Get_dates_daily.py 

# Construct list of dates.
cd $WORKDIR || exit 1
DT=1 # Day
./Get_dates_daily.py $STARTDATE $ENDDATE $DT || exit 1
if [ ! -e "$WORKDIR/Datetime.dat" ] ; then
    echo "ERROR, Datetime.dat not generated!"
    exit 1
fi

# Loop through the dates, fetch date specific files, and run GOCART2WRF.
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

    # Get inst3_3d_aer_Nv file.
    file=$MERRANAME.inst3_3d_aer_Nv.$YEAR$MONTH$DAY.nc4
    if [ ! -e "./$MERRANAME/stage/Y$YEAR/M$MONTH" ] ; then
	mkdir -p "./$MERRANAME/stage/Y$YEAR/M$MONTH" || exit 1
    fi
    if [ ! -e "./$MERRANAME/stage/Y$YEAR/M$MONTH/$file" ] ; then
	WEBSITE=https://goldsmr5.gesdisc.eosdis.nasa.gov
	path=data/MERRA2/M2I3NVAER.5.12.4/$YEAR/$MONTH/$file
	wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies -r -c -nH -nd -np -A nc4  $WEBSITE/$path || exit 1
	mv $file ./$MERRANAME/stage/Y$YEAR/M$MONTH || exit 1
    fi

done

# Now run GOCART2WRF
cd $WORKDIR || exit 1
cat > namelist.gocart2wrf <<EOF
&wrf
  max_dom=$MAXDOM,
  wrf_dir="$WORKDIR",
/

&gocart_shared
  gocart_format=5,
  gocart_dir="$WORKDIR",
  gocart_source="MERRA2",
/
EOF

./gocart2wrf.x || exit 1

# Move output files
if [ ! -e ./wrfbdy_d01 ] ; then
    echo "ERROR, wrfbdy_d01 is not in $WORKDIR!"
    exit 1
fi
mv ./wrfbdy_d01 wrfbdy_d01.gocart2wrf || exit 1
for (( x=1; x <= $MAXDOM; x += 1 )) ; do
    printf -v domain %2.2d $x
    if [ ! -e ./wrfinput_d$domain ] ; then
	echo "ERROR, wrfinput_d$domain is not in $WORKDIR!"
	exit 1
    fi
    mv ./wrfinput_d$domain ./wrfinput_d$domain.gocart2wrf || exit 1
done

# Clean up MERRA2 files
for MERRANAME in MERRA2_100 MERRA2_200 MERRA2_300 MERRA2_400 ; do
    if [ -e "$WORKDIR/$MERRANAME" ] ; then
	rm -rf "$WORKDIR/$MERRANAME"
    fi
done

# The End
echo "GOCART2WRF output files are in $WORKDIR"
echo "Completed MERRA2 GOCART processing"
exit 0
