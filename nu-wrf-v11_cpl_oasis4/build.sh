#!/usr/bin/env bash

export NUWRFDIR=`pwd`

function usage {
   cat $NUWRFDIR/scripts/python/build/README
}

if [ "$#" -eq 0 ]; then
   usage
   exit
fi

source $NUWRFDIR/scripts/other/set_module_env.bash

# Clean up from previous build
rm -f make.log
rm -f .build_failed

cd scripts/python/build
python main.py "$@"
