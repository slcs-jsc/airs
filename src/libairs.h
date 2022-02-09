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

} ret_t;

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

  /*! Fit [K]. */
  double fit[WX][WY];

  /*! Variance [K]. */
  double var[WX][WY];

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

/*! Get day of year from date. */
void day2doy(
  int year,
  int mon,
  int day,
  int *doy);

/*! Get date from day of year. */
void doy2day(
  int year,
  int doy,
  int *mon,
  int *day);

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
  ret_t * ret);

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
  ret_t * ret,
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
