#! /bin/bash

# Get absolute path...
target=$(rm -rf build && mkdir -p build && cd build && pwd)

# Prepare directories...
mkdir -p $target/src $target/bin $target/include $target/lib $target/man/man1 \
    && cp *tar.bz2 $target/src \
    && cd $target/src \
    && for f in *tar.bz2 ; do tar xvjf $f ; done \
	|| exit

# Build only AIRS reader...
if [ "$1" = "-airs" ] ; then
    cd $target/src/airs-v6 \
	&& export CFLAGS="-I$target/include -I/usr/include/hdf -I/usr/include/x86_64-linux-gnu/hdf" \
	&& make \
	&& cp *.h $target/include \
	&& cp lib* $target/lib
    exit
fi

# libjpeg...
dir=jpeg-6b
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make && make check && make install \
    && cp $target/src/$dir/libjpeg.a $target/lib \
    && cp $target/src/$dir/*.h $target/include \
	|| exit

# zlib...
dir=zlib-1.3.1
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
	|| exit

# szip...
dir=szip-2.1.1
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
	|| exit

# tirpc...
dir=libtirpc-1.3.3
cd $target/src/$dir \
    && ./configure --prefix=$target --disable-gssapi \
    && make -j && make check && make install \
	|| exit

# libnsl...
dir=libnsl
cd $target/src/$dir \
    && ./autogen.sh && ./configure --prefix=$target \
    && make -j && make check && make install \
	|| exit

# HDF4...
dir=hdf-4.2.16
cd $target/src/$dir \
    && ./configure --prefix=$target --disable-netcdf --disable-fortran \
		   --disable-hdf4-xdr --with-jpeg=$target \
    && make -j && make check && make install \
	|| exit

# HDF-EOS...
cd $target/src/hdfeos && bin/INSTALL-HDFEOS <<EOF
$target/lib
$target/include
EOF
cp include/* $target/include
cp lib/linux/* $target/lib

# AIRS reader...
cd $target/src/airs-v6 \
    && export CFLAGS=-I$target/include \
    && make \
    && cp *.h $target/include \
    && cp lib* $target/lib

# HDF5...
dir=hdf5-1.14.4-3
cd $target/src/$dir \
    && ./configure --prefix=$target --with-zlib=$target --with-szlib=$target --enable-hl --disable-fortran \
    && make -j && make check && make install && make clean \
	|| exit

# netCDF...
dir=netcdf-c-4.9.2
cd $target/src/$dir \
    && CPPFLAGS=-I$target/include LDFLAGS=-L$target/lib ./configure --prefix=$target --disable-dap --disable-byterange --disable-nczarr --disable-libxml2 \
    && make -j && make install && make clean \
	|| exit

# GSL...
dir=gsl-2.7.1
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
	|| exit
