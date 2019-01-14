#! /bin/bash

# ----------------------------------------------------------------------
function infomsg {
# ----------------------------------------------------------------------
    
    echo
    echo "============================================================"
    echo "Compile: $1"
    echo "============================================================"
    echo
}

# ----------------------------------------------------------------------
# Main...
# ----------------------------------------------------------------------

# Check arguments...
if [ $# -ne 1 ] ; then
    echo "usage: $0 <target-dir>"
    exit
fi

# Set number of cores...
THREADS=$(cat /proc/cpuinfo | grep processor | wc -l)

# Get absolute path...
target=$(mkdir -p $1 && cd $1 && pwd)

# Prepare directories...
mkdir -p $target/src $target/bin $target/lib $target/man/man1 \
    && cp *tar.gz $target/src \
    && cd $target/src \
    && for f in $(ls *tar.gz) ; do tar xvzf $f ; done \
    || exit

# libjpeg...
dir=jpeg-6b
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make && make check && make install \
    && cp $target/src/$dir/libjpeg.a $target/lib \
    || exit

# zlib...
dir=zlib-1.2.3
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j$THREADS && make check && make install \
    || exit

# HDF4...
dir=HDF4.2r4_fix
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target --disable-netcdf --disable-fortran \
    && make -j$THREADS && make check && make install \
    || exit

# HDF-EOS...
infomsg hdfeos
cd $target/src/hdfeos && bin/INSTALL-HDFEOS <<EOF
$target/lib
$target/include
EOF
cp include/* $target/include
cp lib/linux/* $target/lib

# AIRS reader...
infomsg airs-v6
cd $target/src/airs-v6 \
    && export CFLAGS=-I$target/include \
    && make \
    && cp *.h $target/include \
    && cp lib* $target/lib

# netCDF...
dir=netcdf-4.1.2
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target --enable-c-only --disable-dap \
    && make -j$THREADS && make check && make install \
    || exit
