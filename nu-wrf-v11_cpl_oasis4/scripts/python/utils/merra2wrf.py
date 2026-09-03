#!/usr/bin/env python
"""
This script fetches GEOS MERRA2 data and process it for MERRA2WRF.
- Roughly translated from Run_MERRA2.csh script (by Jossy Jacob) -

Try:

  $ merra2wrf.py -h

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
        "-o",
        "--out_dir",
        metavar="out_dir",
        type=str,
        default="./",
        help="Directory where to put the data."
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
        help="Directory where NU-WRF is installed."
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose output")

    return parser.parse_args()


def get_file(merra_server, stage_prefix, base_name, outdir):
    file_path = os.path.join(merra_server, stage_prefix+'/'+base_name)
    utils.copy_file(file_path, outdir)


def main():

    if sys.version_info[0] == 2 and sys.version_info[1] < 7:
        print("Python version " + python_version() + " is not supported.")
        sys.exit()

    # Process options
    args = parse_args()
    start_date = str(args.start_date)
    end_date = str(args.end_date)

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

    print("Start MERRA2WRF conversion...")
    start_time = time.time()

    if not os.path.isfile(nuwrf_dir+'/utils/bin/merra2wrf.x'):
        print(" +++ Error: merra2wrf.x does not exist +++")
        print(" Did you build NU-WRF?")
        sys.exit()

    # The MERRA2WRF application:
    utils.force_symlink(nuwrf_dir +
                        '/utils/bin/merra2wrf.x',
                        './merra2wrf.x')

    print("start_time: "+start_date)
    print("end_time:   "+end_date)
    print("out_dir:    "+out_dir)
    print("nuwrf_dir:  "+nuwrf_dir)

    # Initialize directories
    merra_server = '/discover/nobackup/projects/gmao/merra2/data/products/'

    work_dir = out_dir

    if out_dir == "./":
        work_dir = os.getcwd()

    namelist_out_dir = out_dir + '/data'
    namelist_out_dir1 = out_dir + '/data/merra2wrf'
    outdir = work_dir + '/data'
    outdir1 = out_dir + '/data/merra2wrf'
    utils.mkdir_p(outdir)
    utils.mkdir_p(outdir1)

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
        dd = si.strftime("%Y")+si.strftime("%m")+si.strftime("%d")
        dates.append(dd)
        si = si + datetime.timedelta(days=dt)

    for d in dates:
        year = d[0:4]
        month = d[4:6]
        day = d[6:]
        print('date: ', year, month, day)

        stage_prefix = 'MERRA2_100'
        if 1991 < int(year) < 2001:
            stage_prefix = 'MERRA2_200'
        if 2000 < int(year) < 2011:
            stage_prefix = 'MERRA2_300'
        if int(year) > 2010:
            stage_prefix = 'MERRA2_400'

        print('getting files...')
        os.chdir(outdir)

        file_name1 = stage_prefix+'.const_2d_asm_Nx.00000000.nc4'
        get_file(merra_server, stage_prefix, file_name1, outdir)

        file_name2 = stage_prefix+'.inst6_3d_ana_Nv.'+year+month+day+'.nc4'
        get_file(merra_server,
                 stage_prefix+'/Y'+year+'/M'+month+'/',
                 file_name2, outdir)

        file_name3 = stage_prefix+'.inst6_3d_ana_Np.'+year+month+day+'.nc4'
        get_file(merra_server,
                 stage_prefix+'/Y'+year+'/M'+month+'/',
                 file_name3, outdir)

        file_name4 = stage_prefix+'.tavg1_2d_slv_Nx.'+year+month+day+'.nc4'
        get_file(merra_server,
                 stage_prefix+'/Y'+year+'/M'+month+'/',
                 file_name4, outdir)

        file_name5 = stage_prefix+'.tavg1_2d_ocn_Nx.'+year+month+day+'.nc4'
        get_file(merra_server,
                 stage_prefix+'/Y'+year+'/M'+month+'/',
                 file_name5, outdir)

        file_name6 = stage_prefix+'.tavg1_2d_lnd_Nx.'+year+month+day+'.nc4'
        get_file(merra_server,
                 stage_prefix+'/Y'+year+'/M'+month+'/',
                 file_name6, outdir)

        merradate = year+'-'+month+'-'+day
        os.chdir(work_dir)

        # Create namelist file
        with open('namelist.merra2wrf_daily', 'w') as f:
            foo = \
                "&input\n"\
                "! Directory to write output\n"\
                "outputDirectory = '"+namelist_out_dir1+"',\n"\
                "! Directory with input MERRA files\n"\
                "merraDirectory = '"+namelist_out_dir+"',\n"\
                "! Format and name of const_2d_asm_Nx file\n"\
                "merraFormat_const_2d_asm_Nx = 2,\n"\
                "merraFile_const_2d_asm_Nx = '"+file_name1+"',\n"\
                "! Number of days to process.\n"\
                "numberOfDays = 1,\n"\
                "! Dates of each day being processed (YYYY-MM-DD)\n"\
                "merraDates(1) = '"+merradate+"',\n"\
                "! Format and Names of inst6_3d_ana_Nv files.\n"\
                "merraFormat_inst6_3d_ana_Nv = 2,\n"\
                "merraFiles_inst6_3d_ana_Nv(1) = '"+file_name2+"',\n"\
                "! Names of inst6_3d_ana_Np files.\n"\
                "merraFormat_inst6_3d_ana_Np = 2,\n"\
                "merraFiles_inst6_3d_ana_Np(1) = '"+file_name3+"',\n"\
                "! Names of tavg1_2d_slv_Nx files.\n"\
                "merraFormat_tavg1_2d_slv_Nx = 2,\n"\
                "merraFiles_tavg1_2d_slv_Nx(1) = '"+file_name4+"',\n"\
                "! Names of tavg1_2d_ocn_Nx files.\n"\
                "merraFormat_tavg1_2d_ocn_Nx = 2,\n"\
                "merraFiles_tavg1_2d_ocn_Nx(1) = '"+file_name5+"',\n"\
                "! Names of tavg1_2d_lnd_Nx files.\n"\
                "merraFormat_tavg1_2d_lnd_Nx = 2,\n"\
                "merraFiles_tavg1_2d_lnd_Nx(1) = '"+file_name6+"',\n"\
                "/\n"
            f.write(foo)
        print("running merra2wrf...")
        if verbose:
            utils.sp_call('./merra2wrf.x namelist.merra2wrf_daily')
        else:
            utils.run_shell_command('./merra2wrf.x namelist.merra2wrf_daily')

        # Clean up
        os.remove(outdir+'/'+file_name2)
        os.remove(outdir+'/'+file_name3)
        os.remove(outdir+'/'+file_name4)
        os.remove(outdir+'/'+file_name5)
        os.remove(outdir+'/'+file_name6)

    os.remove(outdir+'/'+file_name1)

    end_time = time.time() - start_time
    print("Time taken = %f" % end_time)
    print(sys.argv[0] + " is DONE.")


if __name__ == "__main__":
    main()
