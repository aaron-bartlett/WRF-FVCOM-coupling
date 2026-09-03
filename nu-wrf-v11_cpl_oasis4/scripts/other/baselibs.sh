#!/usr/bin/env bash
#
# This script installs NU-WRF library dependencies.
# Installation options are set in baselibs.cfg.
#
# Asumptions:
#
# 1) Various paths are hardwired for NCCS's discover system but
#    can be over-ridden by optionally re-defining two environment
#    variables:
#    LIBDIR_TAG: directory where libraries are installed
#    LIBDIR: top-level path for TAG installations
#    URL_LOCAL: location of source code tarballs
#    Examples (these are the defaults for the gnu-openmpi build):
#      LIBDIR_TAG=/discover/nobackup/projects/nu-wrf/lib/gnu-openmpi
#      LIBDIR=/discover/nobackup/projects/nu-wrf/lib
#      URL_LOCAL=/discover/nobackup/projects/nu-wrf/lib/tarballs
#
# 2) Most source code is (can be) obtained from web urls using wget.
#    Some code is stored locally (see URL_LOCAL in (1)).
#
# Note that installation must occur in the order provided in this
# script to resolve inter-library dependencies.

#set -o nounset
#set -o xtrace
#set -o errexit

# Set magic variables for current file & dir
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
__file="${__dir}/$(basename "${BASH_SOURCE[0]}")"
__base="$(basename ${__file} .sh)"
__root="$(cd "$(dirname "${__dir}")" && pwd)" # <-- change this as it depends on your app

arg1="${1:-}"

usage() {
  echo "Usage: $0 [-c|--compiler] option [mp|--mpi] option"
  echo "Possibilities:"
  echo "      $0 -c intel -m sgimpt (intelmpi, mvapich2)"
  echo "      $0 -c gnu -m sgimpt (openmpi,  mvapich2)"
  echo "      $0 -c nag -m mvapich2"
  echo "      $0 -c pgi -m openmpi"
  exit 0
}

MACHINE=$(uname -n)
if [[ "${MACHINE}" =~ "discover" || "${MACHINE}" =~ "borg" ]]; then
   export platform="discover"
else
   export platform="pleiades"
fi

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -h|--help)
    usage
    exit 0
    ;;
    -c|--compiler)
    COMPILER_VENDOR="$2"
    shift # past argument
    ;;
    -m|--mpi)
    MPI_VENDOR="$2"
    shift # past argument
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done
export COMPILER_VENDOR
export MPI_VENDOR

if [ -z "$COMPILER_VENDOR" ] || [ -z "$MPI_VENDOR" ]; then
    usage
fi

source ./set_module_env.bash
module list

mkdir -p $LIBDIR_TAG
# Make sure to revise baselibs.cfg for correct options.
source baselibs.cfg

function get_lib() {
    PACKAGE=$1
    URL=$2
    VER=$3
    SKIP=$4
    
    local PACKAGE_PATH=$LIBDIR_TAG/$PACKAGE
    # If we pass this argument then SKIP removal of previous installation
    if [ -z $SKIP ]; then
        if [ -e $PACKAGE_PATH ]; then
            rm -rf $PACKAGE_PATH || exit 1
        fi
    fi
    cd $LIBDIR

    rm -rf ${PACKAGE}_build
    # This is the build directory:
    mkdir ${PACKAGE}_build

    # Using a stored tarball
    if [[ $URL =~ /discover/nobackup ]]; then
        cp $URL/$VER .
     else
        # Get source code from URL
        wget -q $URL/$VER
        rc=$?
        if [[ $rc != 0 ]]; then
            echo 'ERROR in:'
            echo wget $URL/$VER
            exit $rc;
        fi
    fi
    # Deal with various formats. 
    if [[ $VER =~ \.Z$ ]]; then
      uncompress $VER
      tarball=`echo $VER | sed -e 's/\.[^.]*$//'`
      tar xf $tarball -C  ${PACKAGE}_build --strip-components 2
    elif [[ $VER =~ \.tar$ ]]; then
      tar xf $VER -C ${PACKAGE}_build --strip-components 1
    elif [[ $VER =~ \.xz$ ]]; then
      tar xJf $VER -C ${PACKAGE}_build --strip-components 1
    elif [[ $VER =~ \.gz$ ]]; then
      tar xfz $VER -C ${PACKAGE}_build --strip-components 1
    elif [[ $VER =~ \.zip$ ]]; then
      unzipdir=`echo $VER | sed -e 's/\.[^.]*$//'`
      rm -rf ${PACKAGE}_build
      unzip -q $VER && mv ${unzipdir} ${PACKAGE}_build
    else
      echo "ERROR: Unknown format"; exit 1
    fi
    cd ${PACKAGE}_build
}


function cleanup() {
    local PACKAGE=$1
    cd $LIBDIR
    mkdir -p $LIBDIR/${COMPILER_VENDOR}_logs
    mv *.install *.make *.config  $LIBDIR/${COMPILER_VENDOR}_logs > /dev/null 2>&1
    rm -f *.tar *.gz $VER
    rm -rf ${PACKAGE}_build
}


function check_prereq() {
    dir=$1
    name=${dir##*/}
    prereq_result=0
    if [ ! -d $dir ]; then
        echo " *** Package $name installation failed. ***"
        prereq_result=1
    fi
}


function check_install() {
    dir=$1
    name=${dir##*/}
    install_result=0
    if [ ! -z "$2" ]; then
       lib=$2
    else
       lib=$name
    fi
    if [ ! -d $dir ] || [ ! -f $dir/lib/lib"$lib".a ]; then
       echo " *** Package $name installation failed. ***"
       install_result=1
    else
       echo " --- $name installation was successful ---"
    fi
}

function check_pre_install() {
    dir=$1
    name=${dir##*/}
    if [ ! -z "$2" ]; then
       lib=$2
    else
       lib=$name
    fi
    need_install="no"
    if [[ "$name$" =~ "esmf" ]]; then
       rc=`find $dir -name libesmf.a`
       if [ "$rc" == "" ]; then
         need_install='yes'
       fi
    else
       if [ ! -f $dir/lib/lib"$lib".a ]; then
         need_install="yes"
       fi
    fi
}


if $zlib_flag; then
    PACKAGE=zlib
    echo " --- $PACKAGE ---"

    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX "z"
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $zlib_url $zlib_version

      echo '  configuring...'
      ./configure --64 --static --prefix=$INSTALL_PREFIX > $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1 
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1 
    
      check_install $INSTALL_PREFIX "z"
      cleanup $PACKAGE $zlib_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
    
fi

if $png_flag; then
    PACKAGE=png
    echo " --- $PACKAGE ---"

    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE    
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $png_url $png_version
    
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      ./configure --prefix=$INSTALL_PREFIX --disable-shared >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 >  $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $png_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $jpeg_flag; then
    PACKAGE=jpeg
    echo " --- $PACKAGE ---"

    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE    
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $jpeg_url $jpeg_version

      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      # Version 6b:
      #./configure --prefix=$INSTALL_PREFIX --disable-f90 --disable-shared >> $LIBDIR/$PACKAGE.config 2>&1
      # >6b
      ./configure --prefix=$INSTALL_PREFIX --enable-shared >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      # For version 6b only (originally set by EMK):
      # Make install only handles binaries and man pages. Also run
      # install-lib for libraries and include files.
      #mkdir -p $INSTALL_PREFIX/bin
      #mkdir -p $INSTALL_PREFIX/man/man1
      #mkdir -p $INSTALL_PREFIX/lib 
      #mkdir -p $INSTALL_PREFIX/include 
      #make install > $LIBDIR/$PACKAGE.install 2>&1
      #make install-lib > $LIBDIR/$PACKAGE.install 2>&1
      # For newer versions ( >6b ) use:
      make install > $LIBDIR/$PACKAGE.install 2>&1

      check_install $INSTALL_PREFIX 
      [ $install_result -eq 0 ] && cleanup $PACKAGE $jpeg_version
      
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi


if $hdf4_flag; then
    PACKAGE=hdf4
    echo " --- $PACKAGE ---"

    # Dependencies: (omitting jpeg)
    #JPEG_PATH=$LIBDIR_TAG/jpeg
    #check_prereq $JPEG_PATH

    #ZLIB_PATH=$LIBDIR_TAG/zlib
    #check_prereq $ZLIB_PATH

    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX "df"
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $hdf4_url $hdf4_version
    
      if [ "$COMPILER_VENDOR" == "gnu" ]; then
        cp $PATCHES/hdf4-linux-gnu $LIBDIR/${PACKAGE}_build/config
      fi
      export PROD_FLAGS=""
      export LDFLAGS="-lm" 
      # The HDF4 libraries include an implementation of the netcdf api which can
      # access hdf files. If building with HDF4 and NetCDF it is necessary to build
      # the HDF library with this disabled. So, we need to add "--disable-netcdf"
      # to "configure" and HDF4 will move its embedded NetCDF functions in a
      # different private namespace to avoid name clashes.
      echo '  configuring...'
      #./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      CFLAGS="$CFLAGS" PROD_CFLAGS="$PROD_FLAGS" \
        ./configure --disable-shared \
        --prefix=$INSTALL_PREFIX \
        --enable-fortran \
        --enable-production \
        --disable-netcdf \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  testing...'
      make -j 2 check > $LIBDIR/$PACKAGE.check 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX "df"
      [ $install_result -eq 0 ] && cleanup $PACKAGE $hdf4_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $hdfeos_flag; then
    PACKAGE=hdfeos
    echo " --- $PACKAGE ---"

    # Dependencies: (omitting jpeg)
    #ZLIB_PATH=$LIBDIR_TAG/zlib
    #check_prereq $ZLIB_PATH

    #JPEG_PATH=$LIBDIR_TAG/jpeg
    #check_prereq $JPEG_PATH

    HDF4_PATH=$LIBDIR_TAG/hdf4
    check_prereq $HDF4_PATH
    
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $hdfeos_url $hdfeos_version
      
      export H4CC="$HDF4_PATH/bin/h4cc -Df2cfortran"
      # One of the HDFEOS makefiles is greatly out of date and 
      # hardwired to use gcc and g77 (!!!). We replace it with a
      # patched version.
      cp -f $PATCHES/hdfeos.makefile.utils \
          $LIBDIR/${PACKAGE}_build/util/makefile || exit 1
      export HDFEOS_CC="$HDF4_PATH/bin/h4cc -Df2cFortran"
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      CFLAGS="-Df2cFortran" ./configure --disable-shared \
        --prefix=$INSTALL_PREFIX \
	--enable-install-include \
        --with-hdf4=$HDF4_PATH \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 CC="$HDFEOS_CC" > $LIBDIR/$PACKAGE.make 2>&1
      echo '  testing...'
      make -j 2 CC="$HDFEOS_CC" check > $LIBDIR/$PACKAGE.check 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $hdfeos_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $hdf5_flag; then
    PACKAGE=hdf5
    echo " --- $PACKAGE ---"

    # Dependencies:
    ZLIB_PATH=$LIBDIR_TAG/zlib
    check_prereq $ZLIB_PATH

    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $hdf5_url $hdf5_version
    
      for file in hdf5*; do
	cp $PATCHES/$file $LIBDIR/${PACKAGE}_build/config || exit 1
      done
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      LDFLAGS="-lm" FCFLAGS="$FFLAGS" CFLAGS="$CFLAGS" CXXFLAGS="$CXXFLAGS" \
        ./configure --enable-shared \
        --prefix=$INSTALL_PREFIX \
        --enable-fortran \
        --enable-cxx \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  testing...'
      make -j 2 check > $LIBDIR/$PACKAGE.check 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $hdf5_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $nc4_flag; then
    PACKAGE=netcdf4
    echo " --- $PACKAGE C library ---"
    
    # Dependencies:
    ZLIB_PATH=$LIBDIR_TAG/zlib
    check_prereq $ZLIB_PATH

    HDF5_PATH=$LIBDIR_TAG/hdf5
    check_prereq $HDF5_PATH
   
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX "netcdf"
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $nc4_url $nc4_version

      # Note: --disable-dap gets rid of curl dependency
      export CPPFLAGS="-I$HDF5_PATH/include "
      export LIBS="-L$HDF5_PATH/lib "
      export LD_LIBRARY_PATH=$LIBDIR_TAG/hdf5/lib:$LD_LIBRARY_PATH
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      CPPFLAGS="$CPPFLAGS" LIBS="$LIBS" \
        ./configure --enable-shared \
        --prefix=$INSTALL_PREFIX \
        --disable-dap \
        --enable-static \
        --enable-large-file-tests \
        --enable-netcdf4 \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      #echo '  testing...'
      #make -j 2 check > $LIBDIR/$PACKAGE.check 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      check_install $INSTALL_PREFIX "netcdf"
      [ $install_result -eq 0 ] && cleanup $PACKAGE $nc4_version
      
    else
        echo " --- $PACKAGE C-library is already installed ---"
    fi
fi

if $nc4_flag; then
    PACKAGE=netcdf4
    echo " --- $PACKAGE Fortran library ---"
    
    # Dependencies:
    NC_PATH=$LIBDIR_TAG/netcdf4
    check_prereq $NC_PATH "netcdf"
    
    check_pre_install $INSTALL_PREFIX "netcdff"
    if [ "$need_install" == "yes" ]; then
      export CPPFLAGS="-I$INSTALL_PREFIX/include"
      export CPPFLAGS="$CPPFLAGS -I$HDF5_PATH/include"
      LIB_NETCDF=`$INSTALL_PREFIX/bin/nc-config --libs`
      export LIBS="-L$INSTALL_PREFIX/libs $LIB_NETCDF"
      export LD_LIBRARY_PATH=$LIBDIR_TAG/netcdf4/lib:$LD_LIBRARY_PATH

      echo " --- $PACKAGE fortran library ---"

      get_lib $PACKAGE $nc4fort_url $nc4fort_version "skip"
        
      echo '  configuring...'
      ./configure --help  >> $LIBDIR/$PACKAGE.config 2>&1
      CPPFLAGS="$CPPFLAGS" LIBS="$LIBS" \
        ./configure --enable-shared \
        --prefix=$INSTALL_PREFIX \
        --enable-large-file-tests \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 >> $LIBDIR/$PACKAGE.make 2>&1
      #echo '  testing...'
      #make -j 2 check >> $LIBDIR/$PACKAGE.check 2>&1
      echo '  installing...'
      make install >> $LIBDIR/$PACKAGE.install 2>&1
      check_install $INSTALL_PREFIX "netcdff"
      [ $install_result -eq 0 ] && cleanup $PACKAGE $nc4fort_version
      
    else
        echo " --- $PACKAGE fortran library is already installed ---"
    fi
fi
    
if $nc4_flag; then
    PACKAGE=netcdf4
    echo " --- $PACKAGE C++ library ---"
    
    # Dependencies:
    NC_PATH=$LIBDIR_TAG/netcdf4
    check_prereq $NC_PATH "netcdf"
    
    check_pre_install $INSTALL_PREFIX "netcdf_c++"
    if [ "$need_install" == "yes" ]; then
      export CPPFLAGS="-I$INSTALL_PREFIX/include"
      export CPPFLAGS="$CPPFLAGS -I$HDF5_PATH/include"
      LIB_NETCDF=`$INSTALL_PREFIX/bin/nc-config --libs`
      export LIBS="-L$INSTALL_PREFIX/libs $LIB_NETCDF"
      echo " --- $PACKAGE C++ library ---"
      get_lib $PACKAGE $nc4cxx_url $nc4cxx_version "skip"

      echo '  configuring...'
      ./configure --help  >> $LIBDIR/$PACKAGE.config 2>&1
      CPPFLAGS="$CPPFLAGS" LIBS="$LIBS" \
        ./configure --disable-shared \
        --prefix=$INSTALL_PREFIX \
	--enable-large-file-tests \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 >> $LIBDIR/$PACKAGE.make 2>&1
      #echo '  testing...'
      #make -j 2 check >> $LIBDIR/$PACKAGE.check 2>&1
      echo '  installing...'
      make install >> $LIBDIR/$PACKAGE.install 2>&1
      check_install $INSTALL_PREFIX "netcdf_c++"
      [ $install_result -eq 0 ] && cleanup $PACKAGE $nc4cxx_version
    
      cleanup $PACKAGE $nc4cxx_version   

    else
        echo " --- $PACKAGE C++ library is already installed ---"
    fi
fi

# OTHER LIBRARIES

if $ghostscript_flag; then
    PACKAGE=ghostscript
    echo " --- $PACKAGE ---"

    # Dependencies:
    JPEG_PATH=$LIBDIR_TAG/jpeg
    check_prereq $JPEG_PATH

    ZLIB_PATH=$LIBDIR_TAG/zlib
    check_prereq $ZLIB_PATH

    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE    
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $ghostscript_url $ghostscript_version

      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      ./configure \
        --prefix=$INSTALL_PREFIX \
        --with-jpeg=$JPEG_PATH \
        --with-zlib=$ZLIB_PATH \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1 
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1 
      
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $ghostscript_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $gribapi_flag; then
    PACKAGE=grib_api
    echo " --- $PACKAGE ---"
    
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $gribapi_url $gribapi_version
      
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      ./configure \
        --prefix=$INSTALL_PREFIX \
        --disable-shared \
        --disable-jpeg \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $gribapi_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $freetype_flag; then
    PACKAGE=freetype 
    echo " --- $PACKAGE ---"

    # Dependencies:
    ZLIB_PATH=$LIBDIR_TAG/zlib
    check_prereq $ZLIB_PATH

    PNG_PATH=$LIBDIR_TAG/png
    check_prereq $PNG_PATH

    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $freetype_url $freetype_version
    
      LDFLAGS="-L$PNG_PATH/lib -L$ZLIB_PATH/lib"
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      ./configure \
        --prefix=$INSTALL_PREFIX \
        --with-png=$PNG_PATH \
        --with-zlib=$ZLIB_PATH \
        --with-bzip2=no \
        --enable-shared=no \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1

      # Needed for NU-WRF
      cd $INSTALL_PREFIX/include
      ln -s freetype2 freetype || exit 1
      cd -

      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $freetype_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $jasper_flag; then
    PACKAGE=jasper 
    echo " --- $PACKAGE ---"
    
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $jasper_url $jasper_version 1
    
      echo '  configuring...'
      mkdir foo && cd foo
      #cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX ../ >> $LIBDIR/$PACKAGE.config 2>&1
      cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX -DJAS_ENABLE_SHARED=false ../ >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      cd $INSTALL_PREFIX
      ln -s lib64 lib
      cd -
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $jasper_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $esmf_flag; then
    PACKAGE=esmf
    echo " --- $PACKAGE ---"

    module unload comp/gcc/9.2.0
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $esmf_url $esmf_version
        
      if [ "$COMPILER_VENDOR" == "gnu" ]; then
        if [ "$MPI_VENDOR" == "sgimpt" ]; then
	   cp $PATCHES/build_rules.mk.gfortran.sgimpt.esmf5.2.0rp3 \
              $LIBDIR/${PACKAGE}_build/build_config/Linux.gfortran.default/build_rules.mk \
	   || exit 1
        fi
      else
        if [ "$MPI_VENDOR" == "sgimpt" ]; then
          # Tell SGI MPT underlying compiler name
	  # Use modified build_rules.mk for sgimpt
	  cp $PATCHES/build_rules.mk.intel.sgimpt.esmf6.3.0rp1 \
          $LIBDIR/${PACKAGE}_build/build_config/Linux.intel.default/build_rules.mk \
	    || exit 1
        fi
      fi
      
      export FCFLAGS=$FFLAGS
      export F90FLAGS=$FFLAGS
      
      export ESMF_BOPT="O"
      export ESMF_ABI=64
      export ESMF_DIR=$LIBDIR/${PACKAGE}_build
      export ESMF_INSTALL_PREFIX=$INSTALL_PREFIX

      echo '  configuring/building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      if [ ! -e $INSTALL_PREFIX/lib/libO/Linux.${ESMF_COMPILER}.64.${ESMF_COMM}.default/libesmf.a ]; then
        echo " *** Package $PACKAGE installation failed. ***"
      else
        echo " --- $PACKAGE installation was successful ---"
        cleanup $PACKAGE $esmf_version
      fi
      cleanup $PACKAGE $esmf_version
      rm -f $INSTALL_PREFIX/lib/libO/Linux.${ESMF_COMPILER}.64.${ESMF_COMM}.default/*.so
    else
      echo " --- $PACKAGE is already installed ---"
    fi
    module load comp/gcc/9.2.0
fi

if $gsl_flag; then
    PACKAGE=gsl 
    echo " --- $PACKAGE ---"
 
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $gsl_url $gsl_version
      
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      ./configure \
        --prefix=$INSTALL_PREFIX \
        --disable-shared \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $gsl_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $pixman_flag; then
    PACKAGE=pixman
    echo " --- $PACKAGE ---"

    # Dependencies:
    ZLIB_PATH=$LIBDIR_TAG/zlib
    check_prereq $ZLIB_PATH

    PNG_PATH=$LIBDIR_TAG/png
    check_prereq $PNG_PATH

    INSTALL_PREFIX=$LIBDIR_TAG/pixman
    check_pre_install $INSTALL_PREFIX "pixman-1"
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $pixman_url $pixman_version
    
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      ./configure \
        --prefix=$INSTALL_PREFIX \
        --with-png=$PNG_PATH \
        --with-zlib=$ZLIB_PATH \
        --with-pic \
	--enable-shared=no \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX "pixman-1"
      [ $install_result -eq 0 ] && cleanup $PACKAGE $pixman_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $cairo_flag; then
    PACKAGE=cairo 
    echo " --- $PACKAGE ---"
    
    # Dependencies:
    PIXMAN_PATH=$LIBDIR_TAG/pixman
    check_prereq $PIXMAN_PATH
    
    PNG_PATH=$LIBDIR_TAG/png
    check_prereq $PNG_PATH
    
    ZLIB_PATH=$LIBDIR_TAG/zlib
    check_prereq $ZLIB_PATH
    
    FREETYPE_PATH=$LIBDIR_TAG/freetype
    check_prereq $FREETYPE_PATH
    
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $cairo_url $cairo_version
    
      # Workaround for cairo
      export pixman_CFLAGS="-I$PIXMAN_PATH/include/pixman-1"
      export pixman_LIBS="-L$PIXMAN_PATH/lib -lpixman-1"
      export FREETYPE_LIBS="-L$FREETYPE_PATH/lib -lfreetype"
      export png_LIBS="-L$PNG_PATH/lib -lpng12 -L$ZLIB_PATH/lib -lz"
      CPPFLAGS="-I$PIXMAN_PATH/include/pixman-1"
      CPPFLAGS="$CPPFLAGS -I$FREETYPE_PATH/include/freetype2"
      CPPFLAGS="$CPPFLAGS -I$ZLIB_PATH/include" 
     
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      CPPFLAGS="$CPPFLAGS" \
        ./configure \
        --prefix=$INSTALL_PREFIX \
        --enable-shared=no \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $cairo_version
      unset PKG_CONFIG_PATH
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $bufrlib_flag; then
    PACKAGE=bufrlib
    echo " --- $PACKAGE ---"

    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE    
    check_pre_install $INSTALL_PREFIX "bufr"
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $bufrlib_url $bufrlib_version
      
      echo '  preprocessing...'
      cp $PATCHES/bufrlib-preproc.sh preproc.sh || exit 1
       ./preproc.sh || exit 1

      echo '  building...'
      if [ "$COMPILER_VENDOR" == "gnu" ]; then
         $CC -c -DUNDERSCORE *.c > $LIBDIR/$PACKAGE.make 2>&1
         $FC -c -DUNDERSCORE -fno-second-underscore *.f >> $LIBDIR/$PACKAGE.make 2>&1
      else
         $CC -c -DUNDERSCORE *.c > $LIBDIR/$PACKAGE.make 2>&1
         $FC -c -DUNDERSCORE *.f >> $LIBDIR/$PACKAGE.make 2>&1
      fi
      ar cr libbufr.a *.o
      echo '  installing...'
      mkdir -p $INSTALL_PREFIX/lib $INSTALL_PREFIX/include
      cp -f bufrlib.h $INSTALL_PREFIX/include
      cp -f libbufr.a $INSTALL_PREFIX/lib
      
      check_install $INSTALL_PREFIX "bufr"
      [ $install_result -eq 0 ] && cleanup $PACKAGE $bufrlib_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $g2clib_flag; then
    PACKAGE=g2clib
    echo " --- $PACKAGE ---"

    # Dependencies:
    ZLIB_PATH=$LIBDIR_TAG/zlib
    check_prereq $ZLIB_PATH
    
    PNG_PATH=$LIBDIR_TAG/png
    check_prereq $PNG_PATH
    
    JASPER_PATH=$LIBDIR_TAG/jasper
    check_prereq $JASPER_PATH

    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE    
    check_pre_install $INSTALL_PREFIX "g2c_v1.6.0"
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $g2clib_url $g2clib_version

      INC1=-I$JASPER_PATH/include
      INC2=-I$PNG_PATH/png/include
      INC3=-I$ZLIB_PATH/zlib/include
      export INC_INPUT="$INC1 $INC2 $INC3"
      export CC_INPUT=$CC
      export CFLAGS='-fpic -O2'
      export CFLAGS_INPUT=$CFLAGS

      echo '  building...'
      make 'INC=$(INC_INPUT)' \
          'CC=$(CC_INPUT)' \
          'CFLAGS=$(CFLAGS_INPUT)' \
          >  $LIBDIR/$PACKAGE.make 2>&1 
      echo '  installing...'
      mkdir -p $INSTALL_PREFIX/lib $INSTALL_PREFIX/include
      cp -f *.h $INSTALL_PREFIX/include
      cp -f libg2c_v1.6.0.a $INSTALL_PREFIX/lib
      ln -s $INSTALL_PREFIX/lib/libg2c_v1.6.0.a $INSTALL_PREFIX/lib/libgrib2c.a
    
      check_install $INSTALL_PREFIX "g2c_v1.6.0"
      [ $install_result -eq 0 ] && cleanup $PACKAGE $g2clib_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi


if $ncarg_flag; then
    PACKAGE=ncarg
    echo " --- $PACKAGE ---"
    
    # Dependencies:
    CAIRO_PATH=$LIBDIR_TAG/cairo
    check_prereq $CAIRO_PATH
    
    JPEG_PATH=$LIBDIR_TAG/jpeg
    check_prereq $JPEG_PATH
    
    ZLIB_PATH=$LIBDIR_TAG/zlib
    check_prereq $ZLIB_PATH
    
    HDF4_PATH=$LIBDIR_TAG/hdf4
    check_prereq $HDF4_PATH
    
    NETCDF4_PATH=$LIBDIR_TAG/netcdf4
    check_prereq $NETCDF4_PATH
    
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $ncarg_url $ncarg_version
    
      if [ "$COMPILER_VENDOR" == "gnu" ]; then
        export CFLAGS='-O -ansi -fopenmp -fPIC'
        export FFLAGS='-fPIC -fopenmp'
        export F90FLAGS='-fPIC -fopenmp'
        export CXXFLAGS='-O -ansi -fopenmp -fPIC'
      else
        export CFLAGS='-O -qopenmp -fPIC'
        export FFLAGS='-fPIC -qopenmp'
        export F90FLAGS='-fPIC -qopenmp'
        export CXXFLAGS='-O -qopenmp -fPIC'
      fi

      export CPPFLAGS='-DNDEBUG'
      export NCARG=${PACKAGE}_build
      if [ "$COMPILER_VENDOR" == "gnu" ]; then
	export NCARG_CONFIG=LINUX.64.GNU
      else
	export NCARG_CONFIG=LINUX.64.INTEL
      fi
      cp $PATCHES/ncarg.Configure $LIBDIR/$NCARG/Configure || exit 1
      cp $PATCHES/ncarg.LINUX.64.GNU $LIBDIR/$NCARG/config/LINUX.64.GNU || exit 1
      cd $LIBDIR/$NCARG/config || exit 1
      make -f Makefile.ini  > $LIBDIR/$PACKAGE.config 2>&1
      ./ymake -config `pwd` >> $LIBDIR/$PACKAGE.config 2>&1

      # Pick correct config file
      cp LINUX LINUX.orig || exit 1
      cat $NCARG_CONFIG > LINUX || exit 1
      cd $LIBDIR/$NCARG || exit 1
    
      ./Configure <<EOF
y
n
$INSTALL_PREFIX
y
n
n
/usr/lib64
/usr/include/X11
n
y
EOF
      make Everything  > $LIBDIR/$PACKAGE.make 2>&1
      # No make check exists for NCAR Graphics, so we'll go to install.
      make install  > $LIBDIR/$PACKAGE.install 2>&1

      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $ncarg_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi    
fi

if $boxmg_flag; then
    PACKAGE=boxmg
    echo " --- $PACKAGE ---"
    
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX boxmg_opt_sgl
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $boxmg_url $boxmg_version
      export BOXMG_OS=`uname -s`
      export BOXMG_CPU=`uname -m`
      export BOXMG_CLEVEL=opt
      export BOXMG_COMPILER=intel
      export BOXMG_ARCH=Linux-intel      
      
      echo '  building...'
      make > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      mkdir -p $INSTALL_PREFIX/lib
      mkdir -p $INSTALL_PREFIX/include
      mkdir -p $INSTALL_PREFIX/extras/msg/include
      export BOXMG=$LIBDIR/${PACKAGE}_build
      cp $BOXMG/lib/* $INSTALL_PREFIX/lib
      cp $BOXMG/include/* $INSTALL_PREFIX/include
      cp $BOXMG/extras/msg/include/* $INSTALL_PREFIX/extras/msg/include

      check_install $INSTALL_PREFIX boxmg_opt_sgl
      cleanup $PACKAGE $boxmg_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi

fi

if $gdal_flag; then
    PACKAGE=gdal
    echo " --- $PACKAGE ---"

    JASPER_PATH=$LIBDIR_TAG/jasper
    check_prereq $JASPER_PATH
    
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      #echo "Unload py module"
      #module unload python/GEOSpyD/Ana2019.10_py3.7
      get_lib $PACKAGE $gdal_url $gdal_version

      echo '  configuring...'
      export LDFLAGS="-L$JASPER_PATH/lib64"
      export LIBS="-ljasper"
      export LD_LIBRARY_PATH=$JASPER_PATH/lib64:$LD_LIBRARY_PATH
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      ./configure --prefix=$INSTALL_PREFIX \
                  --with-hide-internal-symbols=yes \
                  --with-gif=internal \
                  --without-webp \
                  --with-libtiff=internal \
                  --with-geotiff=internal \
                  --without-hdf4 \
                  --without-hdf5 \
                  --without-netcdf \
                  --with-jasper=$JASPER_PATH \
                  --with-libjson-c=internal \
                  --without-xml2 >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make -j 2 install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX
      cleanup $PACKAGE $gdal_version
      [ $install_result -eq 0 ] && cleanup $PACKAGE $gdal_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $fortrangis_flag; then
    PACKAGE=fortrangis 
    echo " --- $PACKAGE ---"
    
    # Dependencies:
    GDAL_PATH=$LIBDIR_TAG/gdal
    check_prereq $GDAL_PATH
    
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $fortrangis_url $fortrangis_version
      export CFLAGS="-L$GDAL_PATH/lib"
      export LIBS="-lgdal -lpthread -lstdc++ -ldl"
      export LD_LIBRARY_PATH=$GDAL_PATH/lib:$LD_LIBRARY_PATH

      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      CPPFLAGS="-I$GDAL_PATH/include" \
        LDFLAGS=$LDFLAGS \
        LIBS=$LIBS \
        ./configure \
        --prefix=$INSTALL_PREFIX \
        --disable-shared \
        --enable-gdal \
        --disable-proj \
        --disable-shapelib \
        --disable-doxydoc \
        --disable-readosm \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $fortrangis_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $eccodes_flag; then
    PACKAGE=eccodes 
    echo " --- $PACKAGE ---"
    
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    
    check_pre_install $INSTALL_PREFIX
    if [ "$need_install" == "yes" ]; then
       
      echo '  configuring...'
      mkdir foo && cd foo
      cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX ../ >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      
      check_install $INSTALL_PREFIX
      [ $install_result -eq 0 ] && cleanup $PACKAGE $eccodes_version
    else
      echo " --- $PACKAGE is already installed ---"
    fi
fi

if $nc4_flag_met; then
    PACKAGE=netcdf4_met
    echo " --- $PACKAGE C library ---"
    
    # Dependencies:
    ZLIB_PATH=$LIBDIR_TAG/zlib
    check_prereq $ZLIB_PATH

    HDF5_PATH=$LIBDIR_TAG/hdf5
    check_prereq $HDF5_PATH
   
    INSTALL_PREFIX=$LIBDIR_TAG/$PACKAGE
    check_pre_install $INSTALL_PREFIX "netcdf"
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $nc4_url $nc4_version_met
      export CPPFLAGS="-I$HDF5_PATH/include -I$ZLIB_PATH/include"
      export LIBS="-L$HDF5_PATH/lib -L$ZLIB_PATH/lib -lm"
      export LDFLAGS="$INSTALL_PREFIX/lib"
      echo '  configuring...'
      LDFLAGS="$LDFLAGS" CPPFLAGS="$CPPFLAGS" LIBS="$LIBS" \
        ./configure \
        --prefix=$INSTALL_PREFIX \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 > $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install > $LIBDIR/$PACKAGE.install 2>&1
      check_install $INSTALL_PREFIX "netcdf"
      [ $install_result -eq 0 ] && cleanup $PACKAGE $nc4_version
      
    else
        echo " --- $PACKAGE C-library is already installed ---"
    fi
fi

    
if $nc4cxx_flag_met; then
    PACKAGE=netcdf4_met
    echo " --- $PACKAGE C++ library ---"
    
    # Dependencies:
    NC_PATH=$LIBDIR_TAG/netcdf4
    check_prereq $NC_PATH "netcdf"
    
    check_pre_install $INSTALL_PREFIX "netcdf_c++"
    if [ "$need_install" == "yes" ]; then
      get_lib $PACKAGE $nc4_url $nc4cxx_version_met "skip"
      export CPPFLAGS="$CPPFLAGS -I$HDF5_PATH/include -I$ZLIB_PATH/include"
      export LDFLAGS="$INSTALL_PREFIX/lib"
      echo '  configuring...'
      ./configure --help  > $LIBDIR/$PACKAGE.config 2>&1
      LDFLAGS="$LDFLAGS" CPPFLAGS="$CPPFLAGS" \
        ./configure \
        --prefix=$INSTALL_PREFIX \
        >> $LIBDIR/$PACKAGE.config 2>&1
      echo '  building...'
      make -j 2 >> $LIBDIR/$PACKAGE.make 2>&1
      echo '  installing...'
      make install >> $LIBDIR/$PACKAGE.install 2>&1
      check_install $INSTALL_PREFIX "netcdf_c++"
    
      [ $install_result -eq 0 ] && cleanup $PACKAGE $nc4cxx_version   

    else
        echo " --- $PACKAGE C++ library is already installed ---"
    fi
fi
