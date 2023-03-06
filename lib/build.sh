#! /bin/bash

# Setup...
target=$(rm -rf build && mkdir -p build && cd build && pwd)
threads=$(cat /proc/cpuinfo | grep processor | wc -l)

# Prepare directories...
mkdir -p $target/src $target/bin $target/include $target/lib $target/man/man1 \
    && cp *tar.bz2 $target/src \
    && cd $target/src \
    && for f in $(ls *tar.bz2) ; do tar xvjf $f ; done \
	|| exit

# AIRS reader...
dir=airs-v6
cd $target/src/$dir \
    && cp libairs.* $target/lib \
    && cp airs*.h $target/include

# GSL...
dir=gsl-2.7
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j$threads && make check && make install && make clean \
	|| exit

# HDF5...
dir=hdf5-1.12.1
cd $target/src/$dir \
    && ./configure --prefix=$target --with-zlib=$target --enable-hl \
    && make -j$threads ; make -j$threads && make check && make install && make clean \
	|| exit

# netCDF...
dir=netcdf-c-4.8.1
cd $target/src/$dir \
    && CPPFLAGS=-I$target/include LDFLAGS=-L$target/lib ./configure --prefix=$target --disable-dap --disable-nczarr \
    && make -j$threads && make check && make install && make clean \
	|| exit

# Summary...
echo -e "\n***** gsl-config *****\n"
$target/bin/gsl-config --libs --cflags --version
echo -e "\n***** nc-config *****"
$target/bin/nc-config --all
