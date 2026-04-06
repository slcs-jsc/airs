# Getting Started

## Prerequisites

The project is documented and tested as a Linux build. The standard toolchain is:

- `gcc`
- `make`
- `bash`
- OpenMP support
- an MPI compiler wrapper for `nlte` and `retrieval`

The AIRS programs also depend on:

- GNU Scientific Library (`gsl`)
- netCDF-C
- HDF4
- HDF-EOS
- supporting compression and RPC libraries used by the bundled build

## Clone The Repository

```bash
git clone https://github.com/slcs-jsc/airs.git
cd airs
```

## Build Bundled Libraries

If the required libraries are not already available on the system, build the bundled copies:

```bash
cd libs
./build.sh
```

This populates `libs/build/` with headers and libraries used by the default source build.

## Compile The Tools

```bash
cd src
make
```

By default, `src/Makefile` builds:

- the AIRS analysis utilities listed in `EXC`
- `nlte`
- `retrieval`

The default install destination for `make install` is `../bin`.

## Important Build Flags

- `STATIC=1`: request static linking.
- `OPT=-O3`: optimization level.
- `INFO=1`: emit optimization information.
- `PROF=1`: enable profiling flags.
- `COV=1`: enable coverage instrumentation.

Examples:

```bash
make STATIC=0
make COV=1
make install DESTDIR=../bin
```

## Library Paths

The default include and library search paths in `src/Makefile` assume:

- local bundled libraries in `../libs/build`
- Debian or Ubuntu-style system library paths under `/usr/lib/x86_64-linux-gnu`

If your environment differs, adjust `INCDIR` and `LIBDIR` in `src/Makefile`.

## Build The Documentation

```bash
cd docs
mkdocs build
```

Or through the source `Makefile`:

```bash
cd src
make mkdocs
```

API documentation is generated separately with:

```bash
cd src
make doxygen
```
