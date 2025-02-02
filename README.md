# AIRS Code Collection

The AIRS Code Collection enables data processing and analysis for remote sensing observations captured by NASA's Atmospheric InfraRed Sounder.

[![release (latest by date)](https://img.shields.io/github/v/release/slcs-jsc/airs)](https://github.com/slcs-jsc/airs/releases)
[![commits since latest release (by SemVer)](https://img.shields.io/github/commits-since/slcs-jsc/airs/latest)](https://github.com/slcs-jsc/airs/commits/master)
[![last commit](https://img.shields.io/github/last-commit/slcs-jsc/airs.svg)](https://github.com/slcs-jsc/airs/commits/master)
[![top language](https://img.shields.io/github/languages/top/slcs-jsc/airs.svg)](https://github.com/slcs-jsc/airs/tree/master/src)
[![code size](https://img.shields.io/github/languages/code-size/slcs-jsc/airs.svg)](https://github.com/slcs-jsc/airs/tree/master/src)
[![repo size](https://img.shields.io/github/repo-size/slcs-jsc/airs.svg)](https://github.com/slcs-jsc/airs/tree/master/src)
[![codacy](https://api.codacy.com/project/badge/Grade/a9de7b2239f843b884d2a4eb583726c9)](https://app.codacy.com/gh/slcs-jsc/airs?utm_source=github.com&utm_medium=referral&utm_content=slcs-jsc/airs&utm_campaign=Badge_Grade_Settings)
[![codecov](https://codecov.io/gh/slcs-jsc/airs/branch/master/graph/badge.svg?token=4X6IEHWUBJ)](https://codecov.io/gh/slcs-jsc/airs)
[![tests](https://img.shields.io/github/actions/workflow/status/slcs-jsc/airs/tests.yml?branch=master&label=tests)](https://github.com/slcs-jsc/airs/actions)
[![docs](https://img.shields.io/github/actions/workflow/status/slcs-jsc/airs/docs.yml?branch=master&label=docs)](https://slcs-jsc.github.io/airs)
[![license](https://img.shields.io/github/license/slcs-jsc/airs.svg)](https://github.com/slcs-jsc/airs/blob/master/COPYING)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.14710848.svg)](https://doi.org/10.5281/zenodo.14710848)

## Installation

This documentation describes the installation on a Linux system.
A number of standard tools such as the GNU Compiler Collection (gcc)
and 'make' are required for installation.

Start by downloading the source code from the git repository:

    git clone https://github.com/slcs-jsc/airs.git

Change to the directory airs/ which holds source codes,
libraries, documentation, etc:

    cd airs

The [GNU Scientific Library](https://www.gnu.org/software/gsl) is
required for numerical calculations and the [Unidata netCDF
library](http://www.unidata.ucar.edu/software/netcdf) is needed for
file-I/O.  Furthermore, the HDF4 and HDFEOS libraries are required to
read AIRS data.  Copies of these libraries can be found in the
repository, if they are not available on your system. A script is
provided to build the libraries:

    cd [airs_directory]/libs
    ./build.sh

Next, change to the source directory and edit the Makefile according to
your needs. In particular, check the paths to the libraries
(INCDIR and LIBDIR). Then try to compile the code:

    cd [airs_directory]/src
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

## Further information

The AIRS gravity wave data sets can be found in this repository:

* Hoffmann, Lars, 2021, "AIRS/Aqua Observations of Gravity Waves",
  <https://doi.org/10.26165/JUELICH-DATA/LQAAJA>, Jülich DATA, V1.

These are the main scientific publications that provide information
about the AIRS Code Collection:

* Hoffmann, L., and Alexander, M. J., Retrieval of stratospheric
  temperatures from Atmospheric Infrared Sounder radiance measurements
  for gravity wave studies, J. Geophys. Res., 114, D07105,
  <https://doi.org/10.1029/2008JD011241>, 2009.

* Hoffmann, L., X. Xue, and M. J. Alexander, A global view of
  stratospheric gravity wave hotspots located with Atmospheric
  Infrared Sounder observations, J. Geophys. Res. Atmos., 118,
  416-434, <https://doi.org/10.1029/2012JD018658>, 2013.

* Hoffmann, L., Alexander, M. J., Clerbaux, C., Grimsdell, A. W.,
  Meyer, C. I., Rößler, T., and Tournier, B.: Intercomparison of
  stratospheric gravity wave observations with AIRS and IASI,
  Atmos. Meas. Tech., 7, 4517–4537,
  <https://doi.org/10.5194/amt-7-4517-2014>, 2014.

* Hoffmann, L., Spang, R., Orr, A., Alexander, M. J., Holt, L. A., and
  Stein, O.: A decadal satellite record of gravity wave activity in
  the lower stratosphere to study polar stratospheric cloud formation,
  Atmos. Chem. Phys., 17, 2901-2920,
  <https://doi.org/10.5194/acp-17-2901-2017>, 2017.

More detailed information for users of the AIRS Code Collection is
provided in the [user manual](https://slcs-jsc.github.io/airs).

## License

The AIRS Code Collection is distributed under the GNU GPL v3.
Software libraries distributed along with this software package may have
their own licenses and copyrights, please see corresponding documentation.

## Contact

We are interested in sharing the AIRS Code Collection for research applications.

Please do not hesitate to contact us if you have any further questions:

Dr. Lars Hoffmann

Jülich Supercomputing Centre, Forschungszentrum Jülich  

e-mail: <l.hoffmann@fz-juelich.de>
