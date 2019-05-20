#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of pressure levels for meteorological data. */
#define EP 91

/*! Maximum number of longitudes for meteorological data. */
#define EX 2880

/*! Maximum number of latitudes for meteorological data. */
#define EY 1441

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/*! Scale height [km]. */
#define H0 7.0

/*! Reference pressure [hPa]. */
#define P0 1013.25

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Convert altitude to pressure. */
#define P(z) (P0*exp(-(z)/H0))

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! Meteorological data. */
typedef struct {

  /*! Time [s]. */
  double time;

  /*! Number of longitudes. */
  int nx;

  /*! Number of latitudes. */
  int ny;

  /*! Number of pressure levels. */
  int np;

  /*! Longitude [deg]. */
  double lon[EX];

  /*! Latitude [deg]. */
  double lat[EY];

  /*! Pressure [hPa]. */
  double p[EP];

  /*! Temperature [K]. */
  float t[EX][EY][EP];

} met_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Add variable defintions to netCDF file. */
void addatt(
  int ncid,
  int varid,
  const char *unit,
  const char *long_name);

/*! Auxilary function for interpolation of meteorological data. */
void intpol_met_3d(
  float array[EX][EY][EP],
  int ip,
  int ix,
  int iy,
  double wp,
  double wx,
  double wy,
  double *var);

/*! Spatial interpolation of meteorological data. */
void intpol_met_space(
  met_t * met,
  double p,
  double lon,
  double lat,
  double *t);

/*! Read meteorological data file. */
void read_met(
  char *filename,
  met_t * met);

/*! Extrapolate meteorological data at lower boundary. */
void read_met_extrapolate(
  met_t * met);

/*! Read and convert variable from meteorological data file. */
void read_met_help(
  int ncid,
  char *varname,
  char *varname2,
  met_t * met,
  int np,
  float dest[EX][EY][EP],
  float scl);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  met_t *met;

  static pert_t *pert, *pert2;

  static wave_t wave;

  char pertname[LEN];

  double temp, var_dh, wsum, kp[NSHAPE], kw[NSHAPE];

  int bg_poly_x, itrack, ixtrack, ix, iy, iz, nz,
    ncid, bt_varid, pt_varid, var_varid, dimid[2];

  size_t start[2], count[2];

  /* ------------------------------------------------------------
     Initialize...
     ------------------------------------------------------------ */

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <era.nc> <airs.nc> <kernel.tab>");

  /* Get control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "5", NULL);
  var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "100", NULL);

  /* Alloc... */
  ALLOC(met, met_t, 1);
  ALLOC(pert, pert_t, 1);
  ALLOC(pert2, pert_t, 1);

  /* Read meteorological data... */
  read_met(argv[2], met);

  /* Read AIRS perturbation data... */
  read_pert(argv[3], pertname, pert);

  /* Copy perturbation data... */
  memcpy(pert2, pert, sizeof(pert_t));

  /* Read kernel function... */
  read_shape(argv[4], kp, kw, &nz);
  for (iz = 0; iz < nz; iz++)
    kp[iz] = P(kp[iz]);

  /* ------------------------------------------------------------
     Simulate AIRS data...
     ------------------------------------------------------------ */

  /* Write info... */
  printf("Simulate measurements...\n");

  /* Loop over scans... */
  for (itrack = 0; itrack < pert->ntrack; itrack++) {

    /* Loop over footprints... */
    for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {

      /* Check measured data... */
      if (pert->time[itrack][ixtrack] < 0
	  || pert->lon[itrack][ixtrack] < -180
	  || pert->lon[itrack][ixtrack] > 180
	  || pert->lat[itrack][ixtrack] < -90
	  || pert->lat[itrack][ixtrack] > 90
	  || pert->pt[itrack][ixtrack] < -100
	  || pert->pt[itrack][ixtrack] > 100
	  || !gsl_finite(pert->bt[itrack][ixtrack])
	  || !gsl_finite(pert->pt[itrack][ixtrack])
	  || !gsl_finite(pert->var[itrack][ixtrack])
	  || !gsl_finite(pert->dc[itrack][ixtrack]))
	continue;

      /* Estimate brightness temperature... */
      pert2->bt[itrack][ixtrack] = wsum = 0;
      for (iz = 0; iz < nz; iz++) {
	intpol_met_space(met, kp[iz], pert->lon[itrack][ixtrack],
			 pert->lat[itrack][ixtrack], &temp);
	pert2->bt[itrack][ixtrack] += kw[iz] * temp;
	wsum += kw[iz];
      }
      pert2->bt[itrack][ixtrack] /= wsum;
    }
  }

  /* ------------------------------------------------------------
     Calculate perturbations and variances...
     ------------------------------------------------------------ */

  /* Write info... */
  printf("Get perturbations and variances...\n");

  /* Convert to wave analysis struct... */
  pert2wave(pert2, &wave, 0, pert2->ntrack - 1, 0, pert2->nxtrack - 1);

  /* Estimate background... */
  background_poly(&wave, bg_poly_x, 0);

  /* Compute variance... */
  variance(&wave, var_dh);

  /* Copy data... */
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++) {
      pert2->pt[iy][ix] = wave.pt[ix][iy];
      pert2->var[iy][ix] = wave.var[ix][iy];
    }

  /* ------------------------------------------------------------
     Write to netCDF file...
     ------------------------------------------------------------ */

  /* Write info... */
  printf("Add data to netCDF file...\n");

  /* Open netCDF file... */
  NC(nc_open(argv[3], NC_WRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "NTRACK", &dimid[0]));
  NC(nc_inq_dimid(ncid, "NXTRACK", &dimid[1]));

  /* Enter define mode... */
  NC(nc_redef(ncid));

  /* Add variables... */
  NC(nc_def_var(ncid, "bt_sim", NC_FLOAT, 2, dimid, &bt_varid));
  addatt(ncid, bt_varid, "K", "simulated brightness temperature");
  NC(nc_def_var(ncid, "bt_sim_pt", NC_FLOAT, 2, dimid, &pt_varid));
  addatt(ncid, pt_varid, "K",
	 "simulated brightness temperature perturbation");
  NC(nc_def_var(ncid, "bt_sim_var", NC_FLOAT, 2, dimid, &var_varid));
  addatt(ncid, var_varid, "K^2", "simulated brightness temperature variance");

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Loop over tracks... */
  for (itrack = 0; itrack < pert2->ntrack; itrack++) {

    /* Set array sizes... */
    start[0] = (size_t) itrack;
    start[1] = 0;
    count[0] = 1;
    count[1] = (size_t) pert2->nxtrack;

    /* Write data... */
    NC(nc_put_vara_double(ncid, bt_varid, start, count, pert2->bt[itrack]));
    NC(nc_put_vara_double(ncid, pt_varid, start, count, pert2->pt[itrack]));
    NC(nc_put_vara_double(ncid, var_varid, start, count, pert2->var[itrack]));
  }

  /* Close file... */
  NC(nc_close(ncid));

  /* Free... */
  free(met);
  free(pert);
  free(pert2);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void addatt(
  int ncid,
  int varid,
  const char *unit,
  const char *long_name) {

  /* Set long name... */
  NC(nc_put_att_text(ncid, varid, "long_name", strlen(long_name), long_name));

  /* Set units... */
  NC(nc_put_att_text(ncid, varid, "units", strlen(unit), unit));
}

/*****************************************************************************/

void intpol_met_3d(
  float array[EX][EY][EP],
  int ip,
  int ix,
  int iy,
  double wp,
  double wx,
  double wy,
  double *var) {

  double aux00, aux01, aux10, aux11;

  /* Interpolate vertically... */
  aux00 = wp * (array[ix][iy][ip] - array[ix][iy][ip + 1])
    + array[ix][iy][ip + 1];
  aux01 = wp * (array[ix][iy + 1][ip] - array[ix][iy + 1][ip + 1])
    + array[ix][iy + 1][ip + 1];
  aux10 = wp * (array[ix + 1][iy][ip] - array[ix + 1][iy][ip + 1])
    + array[ix + 1][iy][ip + 1];
  aux11 = wp * (array[ix + 1][iy + 1][ip] - array[ix + 1][iy + 1][ip + 1])
    + array[ix + 1][iy + 1][ip + 1];

  /* Interpolate horizontally... */
  aux00 = wy * (aux00 - aux01) + aux01;
  aux11 = wy * (aux10 - aux11) + aux11;
  *var = wx * (aux00 - aux11) + aux11;
}

/*****************************************************************************/

void intpol_met_space(
  met_t * met,
  double p,
  double lon,
  double lat,
  double *t) {

  double wp, wx, wy;

  int ip, ix, iy;

  /* Check longitude... */
  if (lon < 0)
    lon += 360;

  /* Get indices... */
  ip = locate_irr(met->p, met->np, p);
  ix = locate_reg(met->lon, met->nx, lon);
  iy = locate_reg(met->lat, met->ny, lat);

  /* Get weights... */
  wp = (met->p[ip + 1] - p) / (met->p[ip + 1] - met->p[ip]);
  wx = (met->lon[ix + 1] - lon) / (met->lon[ix + 1] - met->lon[ix]);
  wy = (met->lat[iy + 1] - lat) / (met->lat[iy + 1] - met->lat[iy]);

  /* Interpolate... */
  intpol_met_3d(met->t, ip, ix, iy, wp, wx, wy, t);
}

/*****************************************************************************/

void read_met(
  char *filename,
  met_t * met) {

  int ip, dimid, ncid, varid, year, mon, day, hour;

  size_t np, nx, ny;

  /* Write info... */
  printf("Read meteorological data: %s\n", filename);

  /* Open netCDF file... */
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "lon", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &nx));
  if (nx > EX)
    ERRMSG("Too many longitudes!");

  NC(nc_inq_dimid(ncid, "lat", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &ny));
  if (ny > EY)
    ERRMSG("Too many latitudes!");

  NC(nc_inq_dimid(ncid, "lev", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &np));
  if (np > EP)
    ERRMSG("Too many pressure levels!");

  /* Store dimensions... */
  met->np = (int) np;
  met->nx = (int) nx;
  met->ny = (int) ny;

  /* Read geolocations... */
  NC(nc_inq_varid(ncid, "time", &varid));
  NC(nc_get_var_double(ncid, varid, &met->time));

  NC(nc_inq_varid(ncid, "lev", &varid));
  NC(nc_get_var_double(ncid, varid, met->p));

  NC(nc_inq_varid(ncid, "lon", &varid));
  NC(nc_get_var_double(ncid, varid, met->lon));

  NC(nc_inq_varid(ncid, "lat", &varid));
  NC(nc_get_var_double(ncid, varid, met->lat));

  /* Convert time... */
  year = (int) met->time / 10000;
  met->time -= year * 10000;
  mon = (int) met->time / 100;
  met->time -= mon * 100;
  day = (int) (met->time);
  met->time -= day;
  hour = (int) (met->time * 24.);
  time2jsec(year, mon, day, hour, 0, 0, 0, &met->time);

  /* Check and convert pressure levels... */
  for (ip = 0; ip < met->np; ip++) {
    if (ip > 0 && met->p[ip - 1] > met->p[ip])
      ERRMSG("Pressure levels must be in descending order!");
    met->p[ip] /= 100.;
  }

  /* Read meteorological data... */
  read_met_help(ncid, "T", "t", met, met->np, met->t, 1.0);

  /* Extrapolate data for lower boundary... */
  read_met_extrapolate(met);

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_met_extrapolate(
  met_t * met) {

  int ip, ip0, ix, iy;

  /* Loop over columns... */
  for (ix = 0; ix < met->nx; ix++)
    for (iy = 0; iy < met->ny; iy++) {

      /* Find lowest valid data point... */
      for (ip0 = met->np - 1; ip0 >= 0; ip0--)
	if (!gsl_finite(met->t[ix][iy][ip0]))
	  break;

      /* Extrapolate... */
      for (ip = ip0; ip >= 0; ip--)
	met->t[ix][iy][ip]
	  = (float) LIN(met->p[ip + 1], met->t[ix][iy][ip + 1],
			met->p[ip + 2], met->t[ix][iy][ip + 2], met->p[ip]);
    }
}

/*****************************************************************************/

void read_met_help(
  int ncid,
  char *varname,
  char *varname2,
  met_t * met,
  int np,
  float dest[EX][EY][EP],
  float scl) {

  static float *help;

  int ip, ix, iy, n = 0, varid;

  /* Alloc... */
  ALLOC(help, float,
	EP * EX * EY);

  /* Check if variable exists... */
  if (nc_inq_varid(ncid, varname, &varid) != NC_NOERR)
    if (nc_inq_varid(ncid, varname2, &varid) != NC_NOERR)
      ERRMSG("Cannot read variable!");

  /* Read data... */
  NC(nc_get_var_float(ncid, varid, help));

  /* Copy and check data... */
  for (ip = 0; ip < np; ip++)
    for (iy = 0; iy < met->ny; iy++)
      for (ix = 0; ix < met->nx; ix++) {
	dest[ix][iy][ip] = scl * help[n++];
	if (dest[ix][iy][ip] < -1e10 || dest[ix][iy][ip] > 1e10)
	  dest[ix][iy][ip] = GSL_NAN;
      }

  /* Free... */
  free(help);
}
