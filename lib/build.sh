#! /bin/bash

# Setup...
target=$(rm -rf build && mkdir -p build && cd build && pwd)
threads=$(cat /proc/cpuinfo | grep processor | wc -l)

# Prepare directories...
mkdir -p $target/src $target/bin $target/include $target/lib $target/man/man1 \
    && cp *tar.bz2 $target/src \
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

# tirpc...
dir=libtirpc-1.3.3
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target --disable-gssapi \
    && make -j$THREADS && make check && make install \
    || exit

# HDF4...
dir=hdf-4.2.16
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target --disable-netcdf --disable-fortran \
		   --disable-hdf4-xdr \
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
