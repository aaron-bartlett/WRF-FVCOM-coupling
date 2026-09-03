#!/bin/bash
# Do not change the following settings
#SBATCH -J rip
#SBATCH -t 1:00:00
#SBATCH -o rip.slurm.out
#SBATCH --ntasks=1
# You may edit the following:
#SBATCH -A s0942
##SBATCH --qos=debug
##SBATCH --mail-user=user@nasa.gov

#---------------------- rip script -----------------------------

# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi
# Load config file for modules and paths
source ./common.reg || exit 1

if [ -z "$NUWRFDIR" ] ; then
    echo "ERROR, NUWRFDIR is not defined!"
    exit 1
fi

# Move to work directory
if [ -z "$RUNDIR" ] ; then
    echo "ERROR, RUNDIR is not defined!"
    exit 1
fi

cd $RUNDIR || exit 1

if [ ! -e $NUWRFDIR/scripts/rip/200.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/200.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/200.in $RUNDIR/200.in || exit 1
if [ ! -e $NUWRFDIR/scripts/rip/250.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/250.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/250.in $RUNDIR/250.in || exit 1
if [ ! -e $NUWRFDIR/scripts/rip/300.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/300.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/300.in $RUNDIR/300.in || exit 1
if [ ! -e $NUWRFDIR/scripts/rip/500Vort.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/500Vort.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/500Vort.in $RUNDIR/500Vort.in || exit 1
if [ ! -e $NUWRFDIR/scripts/rip/700RH.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/700RH.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/700RH.in $RUNDIR/700RH.in || exit 1
if [ ! -e $NUWRFDIR/scripts/rip/850.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/850.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/850.in $RUNDIR/850.in || exit 1
if [ ! -e $NUWRFDIR/scripts/rip/sfcThk.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/sfcThk.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/sfcThk.in $RUNDIR/sfcThk.in || exit 1
if [ ! -e $NUWRFDIR/scripts/rip/sfcTUV.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/sfcTUV.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/sfcTUV.in $RUNDIR/sfcTUV.in || exit 1
if [ ! -e $NUWRFDIR/scripts/rip/COMDBZ.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/COMDBZ.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/COMDBZ.in $RUNDIR/COMDBZ.in || exit 1
if [ ! -e $NUWRFDIR/scripts/rip/sfcDBZUV.in ] ; then 
    echo "ERROR, $NUWRFDIR/scripts/rip/sfcDBZUV.in does not exist!"
    exit 1
fi
ln -fs $NUWRFDIR/scripts/rip/sfcDBZUV.in $RUNDIR/sfcDBZUV.in || exit 1
ripfiles="200 250 300 500Vort 700RH 850 sfcThk sfcTUV COMDBZ sfcDBZUV "

# Link rip and rip preprocessor to work directory.
if [ -z "$RIP_ROOT" ] ; then
    echo "ERROR, RIP_ROOT is not defined!"
    exit 1
fi
ln -fs $RIP_ROOT/ripdp_wrfarw ripdp_wrfarw || exit 1
if [ ! -e ripdp_wrfarw ] ; then
    echo "ERROR, ripdp_wrfarw does not exist!"
    exit 1
fi
ln -fs $RIP_ROOT/rip rip || exit 1
if [ ! -e rip ] ; then
    echo "ERROR, rip does not exist!"
    exit 1
fi

# Process each domain
domain_total=0
domain_pass=0
for domain in d01 d02 d03 d04 ; do

    # Count files, and exit for look if no files are found for current domain.
    count=`ls -x -1 | grep wrfout_${domain} | wc -l`

    if [ $count -eq 0 ] ; then
        break
    fi

    let "domain_total = domain_total + 1"

    # Run preprocessor on current domain wrfout files
    files=`ls wrfout_${domain}_*_00:* wrfout_${domain}_*_03:*               wrfout_${domain}_*_06:* wrfout_${domain}_*_09:*               wrfout_${domain}_*_12:* wrfout_${domain}_*_15:*               wrfout_${domain}_*_18:* wrfout_${domain}_*_21:*`
    ./ripdp_wrfarw nuwrf_${domain} all $files || exit 1

    # Now run rip for each rip-execution-name.  Rename the cgm file to
    # prevent overwrites when processing a different domain.
    total=0
    pass=0
    for ripfile in $ripfiles ; do
        ./rip -f nuwrf_${domain} $ripfile || exit 1
        mv ${ripfile}.cgm ${ripfile}_${domain}.cgm || exit 1
        mv ${ripfile}.out ${ripfile}_${domain}.out || exit 1
        test=`grep "We're outta here like Vladimir" ${ripfile}_${domain}.out | wc -l`
        let "total = total + 1"
        if [ $test -eq 1 ] ; then
          let "pass = pass + 1"
        fi
    done
    if [ $total -eq $pass ] ; then
      let "domain_pass = domain_pass + 1"
    fi
    echo "${pass} out of ${total} rip charts succeeded for ${domain} domain" > rip_result_${domain}.out
    rm -fv nuwrf_${domain}_*.00000_* 
done

echo "${domain_pass} out of ${domain_total} domains succeeded" > rip_results.out
if [ $domain_total -eq $domain_pass ] ; then
  echo "Success" >> rip_results.out
fi

# Tidy up logs
mkdir -p rip_logs || exit 1

mv rip_result* rip_logs

# The end
exit 0
