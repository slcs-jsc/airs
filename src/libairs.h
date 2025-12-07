/*
  This file is part of the AIRS Code Collection.
  
  the AIRS Code Collections is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.
  
  The AIRS Code Collection is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with the AIRS Code Collection. If not, see
  <http://www.gnu.org/licenses/>.
  
  Copyright (C) 2019-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  AIRS Code Collection library declarations.
*/

/*! 
  \mainpage

  The AIRS Code Collection enables data processing and analysis
  for remote sensing observations captured by NASA's
  Atmospheric InfraRed Sounder.

  \section Introduction

  The source code of the AIRS Code Collection is available from the
  [git repository](https://github.com/slcs-jsc/airs). Please see the
  [README.md](https://github.com/slcs-jsc/airs/blob/master/README.md)
  in the git repository for introductory information. More information
  can be found in the [user manual](https://slcs-jsc.github.io/airs).
  
  This doxygen manual contains information about the algorithms and
  data structures used in the code. Please refer to the `libairs.h'
  documentation for a first overview.
  
  \section References
  
  For citing the model in scientific publications, please see
  [CITATION.cff](https://github.com/slcs-jsc/airs/blob/master/CITATION.cff).
  
  \section License
  
  The AIRS Code Collection is being develop at the Jülich Supercomputing Centre,
  Forschungszentrum Jülich, Germany.
  
  the AIRS Code Collection is distributed under the terms of the
  [GNU General Public License v3.0](https://github.com/slcs-jsc/airs/blob/master/COPYING).
  
  \section Contributing
  
  We are interested in supporting operational and research
  applications with the AIRS Code Collection.
  
  You can submit bug reports or feature requests on the
  [issue tracker](https://github.com/slcs-jsc/airs/issues).
  
  Proposed code changes and fixes can be submitted as
  [pull requests](https://github.com/slcs-jsc/airs/pulls).
  
  Please do not hesitate to contact us if you have any questions or
  need assistance.
  
  \section Contact
  
  Dr. Lars Hoffmann
  
  Jülich Supercomputing Centre, Forschungszentrum Jülich
  
  e-mail: <l.hoffmann@fz-juelich.de>
*/

#include <netcdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#include <airs_rad_typ.h>
#include <airs_rad_struct.h>
#include <airs_ret_typ.h>
#include <airs_ret_struct.h>
#include "jurassic.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of data sets per granule. */
#define NDS 200000

/*! Maximum number of data points per granule. */
#define NPG 30

/*! Number of AIRS radiance channels (don't change). */
#define L1_NCHAN 34

/*! Along-track size of AIRS radiance granule (don't change). */
#define L1_NTRACK 135

/*! Across-track size of AIRS radiance granule (don't change). */
#define L1_NXTRACK 90

/*! Number of AIRS pressure layers (don't change). */
#define L2_NLAY 27

/*! Along-track size of AIRS retrieval granule (don't change). */
#define L2_NTRACK 45

/*! Across-track size of AIRS retrieval granule (don't change). */
#define L2_NXTRACK 30

/*! Along-track size of perturbation data. */
#define PERT_NTRACK 132000

/*! Across-track size of perturbation data. */
#define PERT_NXTRACK 360

/*! Across-track size of wave analysis data. */
#define WX 300

/*! Along-track size of wave analysis data. */
#define WY 33000

/*! Maximum number of data points for spectral analysis. */
#define PMAX 512

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Execute netCDF library command and check result. */
#define NC(cmd) {				     \
  int nc_result=(cmd);				     \
  if(nc_result!=NC_NOERR)			     \
    ERRMSG("%s", nc_strerror(nc_result));	     \
}

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! AIRS Level-1 data. */
typedef struct {

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[L1_NTRACK][L1_NXTRACK];

  /*! Footprint longitude [deg]. */
  double lon[L1_NTRACK][L1_NXTRACK];

  /*! Footprint latitude [deg]. */
  double lat[L1_NTRACK][L1_NXTRACK];

  /*! Satellite altitude [km]. */
  double sat_z[L1_NTRACK];

  /*! Satellite longitude [deg]. */
  double sat_lon[L1_NTRACK];

  /*! Satellite latitude [deg]. */
  double sat_lat[L1_NTRACK];

  /*! Channel frequencies [cm^-1]. */
  double nu[L1_NCHAN];

  /*! Radiance [W/(m^2 sr cm^-1)]. */
  float rad[L1_NTRACK][L1_NXTRACK][L1_NCHAN];

} airs_l1_t;

/*! AIRS Level-2 data. */
typedef struct {

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[L2_NTRACK][L2_NXTRACK];

  /*! Geopotential height [km]. */
  double z[L2_NTRACK][L2_NXTRACK][L2_NLAY];

  /*! Longitude [deg]. */
  double lon[L2_NTRACK][L2_NXTRACK];

  /*! Latitude [deg]. */
  double lat[L2_NTRACK][L2_NXTRACK];

  /*! Pressure [hPa]. */
  double p[L2_NLAY];

  /*! Temperature [K]. */
  double t[L2_NTRACK][L2_NXTRACK][L2_NLAY];

} airs_l2_t;

/*! Perturbation data. */
typedef struct {

  /*! Number of along-track values. */
  int ntrack;

  /*! Number of across-track values. */
  int nxtrack;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[PERT_NTRACK][PERT_NXTRACK];

  /*! Longitude [deg]. */
  double lon[PERT_NTRACK][PERT_NXTRACK];

  /*! Latitude [deg]. */
  double lat[PERT_NTRACK][PERT_NXTRACK];

  /*! Brightness temperature (8 micron) [K]. */
  double dc[PERT_NTRACK][PERT_NXTRACK];

  /*! Brightness temperature (4 or 15 micron) [K]. */
  double bt[PERT_NTRACK][PERT_NXTRACK];

  /*! Brightness temperature perturbation (4 or 15 micron) [K]. */
  double pt[PERT_NTRACK][PERT_NXTRACK];

  /*! Brightness temperature variance (4 or 15 micron) [K]. */
  double var[PERT_NTRACK][PERT_NXTRACK];

} pert_t;

/*! Retrieval results. */
typedef struct {

  /*! Number of data sets. */
  int nds;

  /*! Number of data points. */
  int np;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[NDS][NPG];

  /*! Altitude [km]. */
  double z[NDS][NPG];

  /*! Longitude [deg]. */
  double lon[NDS][NPG];

  /*! Latitude [deg]. */
  double lat[NDS][NPG];

  /*! Pressure [hPa]. */
  double p[NDS][NPG];

  /*! Temperature [K]. */
  double t[NDS][NPG];

  /*! Temperature (a priori data) [K]. */
  double t_apr[NDS][NPG];

  /*! Temperature (total error) [K]. */
  double t_tot[NDS][NPG];

  /*! Temperature (noise error) [K]. */
  double t_noise[NDS][NPG];

  /*! Temperature (forward model error) [K]. */
  double t_fm[NDS][NPG];

  /*! Temperature (measurement content). */
  double t_cont[NDS][NPG];

  /*! Temperature (resolution). */
  double t_res[NDS][NPG];

  /*! Chi^2. */
  double chisq[NDS];

} retr_t;

/*! Wave analysis data. */
typedef struct {

  /*! Number of across-track values. */
  int nx;

  /*! Number of along-track values. */
  int ny;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time;

  /*! Altitude [km]. */
  double z;

  /*! Longitude [deg]. */
  double lon[WX][WY];

  /*! Latitude [deg]. */
  double lat[WX][WY];

  /*! Across-track distance [km]. */
  double x[WX];

  /*! Along-track distance [km]. */
  double y[WY];

  /*! Temperature [K]. */
  double temp[WX][WY];

  /*! Background [K]. */
  double bg[WX][WY];

  /*! Perturbation [K]. */
  double pt[WX][WY];

  /*! Variance [K]. */
  double var[WX][WY];

  /*! Fit [K]. */
  double fit[WX][WY];

} wave_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Add variable attributes to netCDF file. */
void add_att(
  int ncid,
  int varid,
  const char *unit,
  const char *long_name);

/*! Add variable to netCDF file. */
void add_var(
  int ncid,
  const char *varname,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int *varid,
  int ndims);

/*! Get background based on polynomial fits. */
void background_poly(
  wave_t * wave,
  int dim_x,
  int dim_y);

/*! Get background based on polynomial fits. */
void background_poly_help(
  double *xx,
  double *yy,
  int n,
  int dim);

/*! Smooth background. */
void background_smooth(
  wave_t * wave,
  int npts_x,
  int npts_y);

/*! Set background... */
void create_background(
  wave_t * wave);

/*! Add noise to perturbations and temperatures... */
void create_noise(
  wave_t * wave,
  double nedt);

/*! Add linear wave pattern... */
void create_wave(
  wave_t * wave,
  double amp,
  double lx,
  double ly,
  double phi,
  double fwhm);

/*! Evaluate wave fit... */
void fit_wave(
  wave_t * wave,
  double amp,
  double phi,
  double kx,
  double ky,
  double *chisq);

/*! Calculate 1-D FFT... */
void fft_help(
  double *fcReal,
  double *fcImag,
  int n);

/*! Calculate 2-D FFT... */
void fft(
  wave_t * wave,
  double *Amax,
  double *phimax,
  double *lhmax,
  double *kxmax,
  double *kymax,
  double *alphamax,
  double *betamax,
  char *filename);

/*! Apply Gaussian filter to perturbations... */
void gauss(
  wave_t * wave,
  double fwhm);

/*! Apply Hamming filter to perturbations... */
void hamming(
  wave_t * wave,
  int nit);

/*! Interpolate to regular grid in x-direction. */
void intpol_x(
  wave_t * wave,
  int n);

/*! Apply median filter to perturbations... */
void median(
  wave_t * wave,
  int dx);

/*! Merge wave structs in y-direction. */
void merge_y(
  wave_t * wave1,
  wave_t * wave2);

/*! Estimate noise. */
void noise(
  wave_t * wave,
  double *mu,
  double *sig);

/*! Compute periodogram. */
void period(
  wave_t * wave,
  double lxymax,
  double dlxy,
  double *Amax,
  double *phimax,
  double *lhmax,
  double *kxmax,
  double *kymax,
  double *alphamax,
  double *betamax,
  char *filename);

/*! Convert radiance perturbation data to wave analysis struct. */
void pert2wave(
  pert_t * pert,
  wave_t * wave,
  int track0,
  int track1,
  int xtrack0,
  int xtrack1);

/*! Read AIRS Level-1 data. */
void read_l1(
  char *filename,
  airs_l1_t * l1);

/*! Read AIRS Level-2 data. */
void read_l2(
  char *filename,
  airs_l2_t * l2);

/*! Read radiance perturbation data. */
void read_pert(
  char *filename,
  char *pertname,
  pert_t * pert);

/*! Read AIRS retrieval data. */
void read_retr(
  char *filename,
  retr_t * ret);

/*! Convert array. */
void read_retr_help(
  double *help,
  int nds,
  int np,
  double mat[NDS][NPG]);

/*! Read wave analysis data. */
void read_wave(
  char *filename,
  wave_t * wave);

/*! Convert AIRS radiance data to wave analysis struct. */
void rad2wave(
  airs_rad_gran_t * airs_rad_gran,
  double *nu,
  int nd,
  wave_t * wave);

/*! Convert AIRS retrieval results to wave analysis struct. */
void ret2wave(
  retr_t * ret,
  wave_t * wave,
  int dataset,
  int ip);

/*! Compute local variance. */
void variance(
  wave_t * wave,
  double dh);

/*! Write AIRS Level-1 data. */
void write_l1(
  char *filename,
  airs_l1_t * l1);

/*! Write AIRS Level-2 data. */
void write_l2(
  char *filename,
  airs_l2_t * l2);

/*! Write wave analysis data. */
void write_wave(
  char *filename,
  wave_t * wave);
