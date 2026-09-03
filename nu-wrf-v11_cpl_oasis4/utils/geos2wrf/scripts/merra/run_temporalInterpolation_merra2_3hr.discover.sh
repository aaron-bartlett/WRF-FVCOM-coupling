#!/bin/sh
#SBATCH --job-name=ti
#SBATCH --time=0:10:00
#SBATCH --account s0942
#SBATCH --output ti.slurm.out
#Adjust node, core, and hardware constraints here.
#SBATCH --ntasks=1 --constraint=hasw
#Substitute your e-mail here
##SBATCH --mail-user=user@nasa.gov
##SBATCH --mail-type=ALL
#Set quality of service, if needed.
##SBATCH --qos=debug
#------------------------------------------------------------------------------
# NASA/GSFC, Computational and Information Science and Technology Office,
# Code 606
#------------------------------------------------------------------------------
#                                                                              
# SCRIPT:  run_temporalInterpolation_merra2_3hr.discover.sh
#                                                                              
# AUTHOR:                                                                      
# Eric Kemp, NASA CISTO/SSAI
#
# DESCRIPTION:
# Sample script for running GEOS2WRF temporalInterpolation tool to interpolate
# 6-hourly MERRA2 data to every three hours. The MERRA2 data must first have
# been processed by MERRA2WRF (i.e., temporalInterpolation will operate on
# WPS intermediate format files).
#
#------------------------------------------------------------------------------

# Strings must be YYYYMMDDHH
analStartTimeString=2007011900 # First MERRA2 reanalysis
analEndTimeString=2007012018   # Last MERRA2 reanalysis
analDeltaHour=6                # Hours between MERRA2 analyses

interpStartTimeString=2007011903 # First interpolation date/time
interpEndTimeString=2007012015   # Last interpolation date/time
interpDeltaHour=3                # Hours between interpolation date/time

prefix="MERRA"

inputDir=/discover/nobackup/projects/nu-wrf/members/emkemp/cases_v7/R_GEOS2WRF_TINTERP/input
outputDir=/discover/nobackup/projects/nu-wrf/members/emkemp/cases_v7/R_GEOS2WRF_TINTERP/output

# When a batch script is started, it starts in the user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Load config file for modules and paths
source ./config.discover.sh || exit 1

# Move to work directory.
if [ -z "$WORKDIR" ] ; then
    echo "ERROR, WORKDIR is not defined!"
    exit 1
fi
cd $WORKDIR || exit 1

# Create links to temporalInterpolation
if [ -z "$NUWRFDIR" ] ; then
    echo "ERROR, NUWRFDIR is not defined!"
    exit 1
fi
for exec in temporalInterpolation ; do
    file="$NUWRFDIR/utils/geos2wrf/$exec"
    if [ ! -e "$file" ] ; then
        echo "ERROR, $file does not exist!"
        exit 1
    fi
    ln -fs $file $WORKDIR/$exec || exit 1
done

# Create link to updateDateTime.py
file="$NUWRFDIR/utils/geos2wrf/scripts/util/updateDateTime.py"
if [ ! -e "$file" ] ; then
    echo "ERROR, $file does not exist!"
    exit 1
fi
ln -fs $file $WORKDIR/updateDateTime.py || exit 1

# Loop through each interpolation time.
interpCurrentTimeString=${interpStartTimeString}
while [ "$interpCurrentTimeString" -le "$interpEndTimeString" ] ; do

    icyear=`echo $interpCurrentTimeString | cut -c1-4`
    icmonth=`echo $interpCurrentTimeString | cut -c5-6`
    icday=`echo $interpCurrentTimeString | cut -c7-8`
    ichour=`echo $interpCurrentTimeString | cut -c9-10`
    iTimeStamp="${icyear}-${icmonth}-${icday}_${ichour}"

    # Loop through available MERRA2 date/times and find nearest before and
    # after interpolation. If an exact match is found to the interpolation
    # date/time, no interpolation will be needed.

    if [ "$analStartTimeString" -gt "$interpCurrentTimeString" ] ; then
	echo "ERROR, MERRA2 reanalyses start after interpolation time $interpCurentTimeString!"
	exit 1
    fi
    if [ "$analEndTimeString" -lt "$interpCurrentTimeString" ] ; then
	echo "ERROR, MERRA2 reanalyses end before interpolation time $interpCurrentTimeString!"
	exit 1
    fi

    foundWindow=0
    analCurrentTimeString=${analStartTimeString}
    while [ "$analCurrentTimeString" -lt "$analEndTimeString" ] ; do
	
	analBeforeTimeString=${analCurrentTimeString}
	analAfterTimeString=`./updateDateTime.py $analCurrentTimeString $analDeltaHour`

	abyear=`echo $analBeforeTimeString | cut -c1-4`
	abmonth=`echo $analBeforeTimeString | cut -c5-6`
	abday=`echo $analBeforeTimeString | cut -c7-8`
	abhour=`echo $analBeforeTimeString | cut -c9-10`
	valid=0
	for hour in "00" "06" "12" "18" ; do
	    if [ "$abhour" -eq "$hour" ] ; then
		valid=1
	    fi
	done
	if [ ! "$valid" -eq 1 ] ; then
	    echo "ERROR, invalid hour for MERRA2 data!"
	    echo "abhour = $abhour"
	    echo "Must be 00, 06, 12, or 18!"
	    exit 1
	fi
	
	aayear=`echo $analAfterTimeString | cut -c1-4`
	aamonth=`echo $analAfterTimeString | cut -c5-6`
	aaday=`echo $analAfterTimeString | cut -c7-8`
	aahour=`echo $analAfterTimeString | cut -c9-10`
	valid=0
	for hour in "00" "06" "12" "18" ; do
	    if [ "$aahour" -eq "$hour" ] ; then
		valid=1
	    fi
	done
	if [ ! "$valid" -eq 1 ] ; then
	    echo "ERROR, invalid hour for MERRA2 data!"
	    echo "aahour = $aahour"
	    echo "Must be 00, 06, 12, or 18!"
	    exit 1
	fi
	
	# See if we can skip interpolation.
	if [ "$analBeforeTimeString" -eq "$interpCurrentTimeString" ] ; then
	    echo "No interpolation needed, MERRA2 analysis available at $interpCurrentTimeString"
	    foundWindow=1
	    break
	fi
	if [ "$analAfterTimeString" -eq "$interpCurrentTimeString" ] ; then
	    echo "No interpolation needed, MERRA2 analysis available at $interpCurrentTimeString"
	    foundWindow=1
	    break
	fi

	if [ "$analBeforeTimeString" -lt "$interpCurrentTimeString" ] ; then
	    if [ "$analAfterTimeString" -gt "$interpCurrentTimeString" ] ; then
		foundWindow=1

		echo "Matched for interpolation!"
		echo "Interpolation time: $interpCurrentTimeString"
		echo "First MERRA2 analysis time: $analBeforeTimeString"
		echo "Second MERRA2 analysis time: $analAfterTimeString"

		fileList=""
		for var in PRESSURE HGT TT UU VV RH PSFC PMSL SKINTEMP \
                           SEAICE SOILHGT LANDSEA ; do
		    cat > namelist.temporalInterpolation<<EOF

&all
  fieldName='${var}',
/ 

&input1
  directory1='${inputDir}',
  prefix1='${prefix}',
  year1=${abyear},
  month1=${abmonth},
  day1=${abday},
  hour1=${abhour},
/

&input2
  directory2='${inputDir}',
  prefix2='${prefix}',
  year2=${aayear},
  month2=${aamonth},
  day2=${aaday},
  hour2=${aahour},
/

&output
  directoryOutput='${outputDir}',
  yearOutput=${icyear},
  monthOutput=${icmonth},
  dayOutput=${icday},
  hourOutput=${ichour},
/

EOF
		    # Run temporalInterpolation. Serial execution.
		    ./temporalInterpolation || exit 1

		    fileList="$fileList $outputDir/$var:$iTimeStamp"

		done

		# Consolidate output slab files
		cat $fileList > $outputDir/$prefix:$iTimeStamp || exit 1
		rm $fileList || exit 1

	    fi
	fi

        #----------------------------------------------------------------------
        # Go to next analysis time
        #----------------------------------------------------------------------
	analCurrentTimeString=${analAfterTimeString}

    done

    if [ "$foundWindow" -eq 0 ] ; then
	echo "ERROR, could not interpolate for $interpCurrentTimeString!"
	exit 1
    fi

    #--------------------------------------------------------------------------
    # Go to next interpolation time
    #--------------------------------------------------------------------------
    interpCurrentTimeString=`./updateDateTime.py $interpCurrentTimeString $interpDeltaHour`
done

# The end
echo "Processing completed."
exit 0
