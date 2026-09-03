#!/usr/bin/env python
"""
This script fetches SST data and process it for SST2WRF.
- Roughly translated from Run_SST.csh script (by Jossy Jacob) -

SST data from the following websites will be downloaded and processed.

   http://data.remss.com/sst

Notes for version-4 (as of 9/2020) 
(1) mw_ir (Microwave IR) is available up to 2017
(2) mw (Microwave-only) is available up to 2018

Notes for version-5 (as of 9/2020) 
(1) mw_ir (Microwave IR) is available mid-2002-present
(2) mw (Microwave-only) is available 1998-present

Try:

  $ sst2wrf.py -h

"""

import os
import sys
import time
import datetime
from platform import python_version
import argparse
import shared.utils as utils


def parse_args():
    """Parse command line arguments.
    @returns: populated namespace containing parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("start_date", help="Start date (YYMMDD).")
    parser.add_argument("end_date", help="End date (YYMMDD).")
    parser.add_argument(
        "-v",
        "--version",
        metavar="data_version",
        type=str,
        default="v5",
        help="Data version, either v4 or v5 (default v5)"
    )
    parser.add_argument(
        "-i",
        "--inst_type",
        metavar="inst_type",
        type=str,
        default="mw_ir",
        help="mw_ir or mw (default mw_ir)"
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        metavar="out_dir",
        type=str,
        default="./",
        help="Directory where to put the data (default ./)"
    )
    nuwrf_dir = None
    if os.environ.get("NUWRFDIR") is not None:
        nuwrf_dir = os.environ.get("NUWRFDIR")
    parser.add_argument(
        "-n",
        "--nuwrf_dir",
        metavar="nuwrf_dir",
        type=str,
        default=nuwrf_dir,
        help="Directory where NU-WRF is installed (default in $NUWRFDIR)"
    )
    parser.add_argument("-d", "--verbose", action="store_true", help="verbose output")

    return parser.parse_args()


def main():

    if sys.version_info[0] == 2 and sys.version_info[1] < 7:
        print("Python version " + python_version() + " is not supported.")
        sys.exit()

    # Process options
    args = parse_args()
    start_date = str(args.start_date)
    end_date = str(args.end_date)

    if args.version:
        data_version = str(args.version)

    if args.inst_type:
        inst_type = str(args.inst_type)

    if args.out_dir:
        out_dir = str(args.out_dir)

    if args.nuwrf_dir is None:
        print("Please specify nuwrf_dir")
        sys.exit()
    else:
        nuwrf_dir = str(args.nuwrf_dir)

    verbose = False
    if args.verbose:
        verbose = True

    print("Start SST2WRF conversion...")
    start_time = time.time()

    if not os.path.isfile(nuwrf_dir+'/utils/bin/sst2wrf.x'):
        print(" +++ Error: sst2wrf.x does not exist +++")
        print(" Did you build NU-WRF?")
        sys.exit()

    # The SST2WRF application:
    utils.force_symlink(nuwrf_dir +
                        '/utils/bin/sst2wrf.x',
                        './sst2wrf.x')

    print("start_time: "+start_date)
    print("end_time:   "+end_date)
    print("inst_type:  "+inst_type)
    print("data_vers:  "+data_version)
    print("out_dir:    "+out_dir)
    print("nuwrf_dir:  "+nuwrf_dir)

    # Initialize directories
    data_server = 'http://data.remss.com/sst'

    work_dir = out_dir

    if out_dir == "./":
        work_dir = os.getcwd()

    outdir = work_dir + '/data'
    outdir1 = out_dir + '/data/'+inst_type
    utils.mkdir_p(outdir)
    utils.mkdir_p(outdir1)

    syear = start_date[0:4]
    if int(syear) > 2017 and inst_type == 'mw_ir':
        print('SST2WRF cannot process mw_ir data for years >= ', syear)
        sys.exit()

    # Make a list of dates to process
    dates = list()
    sdate = datetime.date(int(start_date[0:4]),
                          int(start_date[4:6]),
                          int(start_date[6:]))
    edate = datetime.date(int(end_date[0:4]),
                          int(end_date[4:6]),
                          int(end_date[6:]))
    si = sdate
    dt = 1  # daily interval
    while si <= edate:
        tt = si.timetuple()
        cday = tt.tm_yday
        dd = "%3.3d%3.3d\n" % (tt.tm_year, tt.tm_yday)
        dates.append(dd)
        si = si + datetime.timedelta(days=dt)

    year = syear
    prefix = inst_type+'.fusion.'
    if "4" in data_version:
        ftpname = data_server+'/daily_v04.0/'+inst_type+'/'+year+'/'
    else:
        ftpname = data_server+'/daily/'+inst_type+'/v05.0/bmaps/'+year+'/'

    if "4" in data_version:
        if inst_type == 'mw':
            ver = 'rt'
        else:
            ver = 'v04.0'
    else:
        ver = 'v05.0'
    suffix = '.'+ver+'.gz'

    for d in dates:
        day = d[4:7]
        if len(day) == 1:
            day = '00'+day
        elif len(day) == 2:
            day = '0'+day
        print(f'getting files for day #', day)
        os.chdir(outdir)
        fname = prefix+year+'.'+day+suffix
        fpname = ftpname+fname
        if verbose:
            utils.sp_call('wget '+fpname)
            utils.sp_call('gunzip '+fname)
        else:
            utils.run_shell_command('wget '+fpname)
            utils.run_shell_command('gunzip '+fname)
            
        os.chdir(work_dir)
        with open('namelist.sst2wrf', 'w') as f:
            foo =\
                "&input\n"\
                "  instrument = '"+inst_type+"',\n"\
                "  year = "+year+",\n"\
                "  dayOfYear = "+day+",\n"\
                "  version = '"+ver+"',\n"\
                "  inputDirectory = 'data/',\n"\
                "/\n"\
                "&output\n"\
                "  outputDirectory = 'data/"+inst_type+"',\n"\
                "  prefixWPS = 'SSTRSS'\n"\
                "/\n"\
                "&fakeoutput\n"\
                "  numFakeHours = 4,\n"\
                "  fakeHours = 0, 6, 12, 18,\n"\
                "/\n"
            f.write(foo)
        print("running sst2wrf...")
        if verbose:
            utils.sp_call('./sst2wrf.x')
        else:
            utils.run_shell_command('./sst2wrf.x')

    end_time = time.time() - start_time
    print("Time taken = %f" % end_time)
    print(sys.argv[0] + " is DONE.")


if __name__ == "__main__":
    main()
