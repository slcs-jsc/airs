# AIRS Code Collection

This repository provides a collection of codes for the analysis of
observations of NASA's Atmospheric InfraRed Sounder (AIRS).

## Installation

This documentation describes the installation on a Linux system.
A number of standard tools such as the GNU Compiler Collection (gcc)
and 'make' are required for installation.

Start by downloading the source code from the github repository:

    git clone https://github.com/slcs-jsc/airs

Change to the directory airs/ which holds source codes,
libraries, documentation, etc:

    cd airs

The GNU Scientific Library (https://www.gnu.org/software/gsl)
is required for numerical calculations and the Unidata netCDF library
(http://www.unidata.ucar.edu/software/netcdf) is needed for file-I/O.
Furthermore, the HDF4 and HDFEOS libraries are required to read AIRS data.
Copies of these libraries can be found in the repository, if they are
not available on your system. A script is provided to build the libraries:

    cd lib
    ./build.sh

Next, change to the source directory and edit the Makefile according to
your needs. In particular, check the paths to the libraries
(INCDIR and LIBDIR). Then try to compile the code:

    cd ../src
    emacs Makefile
    make

The binaries will be linked statically, i.e., they can be copied to other
machines. Sometimes static compilations causes problems, in particular in
combination with MPI. In this case remove the '-static' flag from the
CFLAGS in the Makefile and compile again.

By default we use rather strict compiler warnings.
All warning messages will be turned into errors and no binaries will be
produced. This behavior is enforced by the flag '-Werror'.

The binaries will remain in the src/ directory.

## Contact

We are interested in sharing the AIRS code for research applications.

Please do not hesitate to contact us if you have any further questions:

Dr. Lars Hoffmann  
Forschungszentrum Jülich  
Jülich Supercomputing Centre  
52425 Jülich  
Germany  

e-mail: l.hoffmann@fz-juelich.de

## License

The AIRS Code Collection is distributed under the GNU GPL v3.
Software libraries distributed along with this software package may have
their own licenses and copyrights, please see corresponding documentation.
