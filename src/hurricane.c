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
  Analyse gravity wave data for tropical cyclones.
*/

#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of storms. */
#define NSTORM 9000

/* Maximum number of observation times. */
#define NTIME 140

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Get storm position at given time... */
int get_storm_pos(
  int nobs,
  double time_wmo[NTIME],
  double lon_wmo[NTIME],
  double lat_wmo[NTIME],
  double wind_wmo[NTIME],
  double pres_wmo[NTIME],
  double t,
  int dt,
  int st,
  double x[3],
  double *wind,
  double *dwind,
  double *pres,
  double *dpres);

/* Read variable from netCDF file... */
void read_var(
  int ncid,
  const char varname[],
  size_t nstorm,
  int nobs[NSTORM],
  double x[NSTORM][NTIME]);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  static FILE *out;

  static char filter[LEN], pertname[LEN], set[LEN];

  static double bt4_mean, bt4_var, bt8_min, dpres, dpresbest, dwind,
    dwindbest, lat_wmo[NSTORM][NTIME], latbest, lon_wmo[NSTORM][NTIME],
    lonbest, lonsat, lonstorm, nedt, pres_wmo[NSTORM][NTIME],
    pres, presbest, r2best = 1e100, wind_wmo[NSTORM][NTIME], wind,
    windbest, time_max_pres[NSTORM], time_max_wind[NSTORM],
    time_wmo[NSTORM][NTIME], timebest, xf[PERT_NTRACK][PERT_NXTRACK][3],
    xs[3], z;

  static int dimid, n, ncid, nobs[NSTORM], varid;

  static size_t nstorm, ntime;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <hurr.tab> <ibtracs.nc>"
	   " <pert1.nc> [<pert2.nc> ...]");

  /* Get control parameters... */
  scan_ctl(argc, argv, "SET", -1, "full", set);
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  scan_ctl(argc, argv, "FILTER", -1, "both", filter);
  const double dt230 = scan_ctl(argc, argv, "DT230", -1, "0.16", NULL);
  const double nu = scan_ctl(argc, argv, "NU", -1, "2345.0", NULL);
  const double rmax = scan_ctl(argc, argv, "RMAX", -1, "500", NULL);
  const int dt = (int) scan_ctl(argc, argv, "DT", -1, "0", NULL);
  const int st = (int) scan_ctl(argc, argv, "ST", -1, "0", NULL);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* ------------------------------------------------------------
     Read hurricane tracks...
     ------------------------------------------------------------ */

  /* Write info... */
  printf("Read hurricane tracks: %s\n", argv[3]);

  /* Open netCDF file... */
  NC(nc_open(argv[3], NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "storm", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &nstorm));
  NC(nc_inq_dimid(ncid, "time", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &ntime));
  if (nstorm > NSTORM)
    ERRMSG("Too many storms!");
  if (ntime > NTIME)
    ERRMSG("Too many time steps!");

  /* Read number of observations per storm... */
  NC(nc_inq_varid(ncid, "numObs", &varid));
  NC(nc_get_var_int(ncid, varid, nobs));

  /* Read data... */
  read_var(ncid, "lat_wmo", nstorm, nobs, lat_wmo);
  read_var(ncid, "lon_wmo", nstorm, nobs, lon_wmo);
  read_var(ncid, "time_wmo", nstorm, nobs, time_wmo);
  read_var(ncid, "wind_wmo", nstorm, nobs, wind_wmo);
  read_var(ncid, "pres_wmo", nstorm, nobs, pres_wmo);

  /* Convert units.. */
  for (size_t istorm = 0; istorm < nstorm; istorm++)
    for (int iobs = 0; iobs < nobs[istorm]; iobs++) {
      time_wmo[istorm][iobs] *= 86400.;
      time_wmo[istorm][iobs] -= 4453401600.00;
      lon_wmo[istorm][iobs] *= 0.01;
      lat_wmo[istorm][iobs] *= 0.01;
      wind_wmo[istorm][iobs] *= 0.0514444;
      pres_wmo[istorm][iobs] *= 0.1;
    }

  /* Check data... */
  for (size_t istorm = 0; istorm < nstorm; istorm++)
    for (int iobs = 0; iobs < nobs[istorm]; iobs++) {
      if (pres_wmo[istorm][iobs] <= 800 || pres_wmo[istorm][iobs] >= 1200)
	pres_wmo[istorm][iobs] = GSL_NAN;
      if (wind_wmo[istorm][iobs] <= 0.1)
	wind_wmo[istorm][iobs] = GSL_NAN;
    }

  /* Find time of maximum intensity (lowest pressure)... */
  for (size_t istorm = 0; istorm < nstorm; istorm++) {
    double pmin = 1e100;
    time_max_pres[istorm] = GSL_NAN;
    for (int iobs = 0; iobs < nobs[istorm]; iobs++)
      if (gsl_finite(pres_wmo[istorm][iobs]) && pres_wmo[istorm][iobs] < pmin) {
	pmin = pres_wmo[istorm][iobs];
	time_max_pres[istorm] = time_wmo[istorm][iobs];
      }
  }

  /* Find time of maximum intensity (maximum wind)... */
  for (size_t istorm = 0; istorm < nstorm; istorm++) {
    double wmax = -1e100;
    time_max_wind[istorm] = GSL_NAN;
    for (int iobs = 0; iobs < nobs[istorm]; iobs++)
      if (gsl_finite(wind_wmo[istorm][iobs]) && wind_wmo[istorm][iobs] > wmax) {
	wmax = wind_wmo[istorm][iobs];
	time_max_wind[istorm] = time_wmo[istorm][iobs];
      }
  }

  /* Close netCDF file... */
  NC(nc_close(ncid));

  /* ------------------------------------------------------------
     Analyze AIRS data...
     ------------------------------------------------------------ */

  /* Create file... */
  printf("Write hurricane data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = storm number\n"
	  "# $2  = storm time since first report [hr]\n"
	  "# $3  = storm time since wind maximum [hr]\n"
	  "# $4  = storm time since pressure minimum [hr]\n"
	  "# $5  = match time [s]\n"
	  "# $6  = match longitude [deg]\n"
	  "# $7  = match latitude [deg]\n"
	  "# $8  = match distance [km]\n"
	  "# $9  = wind speed [m/s]\n"
	  "# $10 = wind speed change [m/s/hr]\n");
  fprintf(out,
	  "# $11 = pressure [hPa]\n"
	  "# $12 = pressure change [hPa/hr]\n"
	  "# $13 = 8.1 micron BT minimum [K]\n"
	  "# $14 = 4.3 micron BT variance [K^2]\n"
	  "# $15 = 4.3 micron BT variance (noise-corrected) [K^2]\n"
	  "# $16 = number of footprints\n\n");

  /* Loop over perturbation files... */
  for (int iarg = 4; iarg < argc; iarg++) {

    /* Read perturbation data... */
    FILE *in;
    if (!(in = fopen(argv[iarg], "r")))
      continue;
    else {
      fclose(in);
      read_pert(argv[iarg], pertname, pert);
    }

    /* Get Cartesian coordinates... */
    for (int itrack2 = 0; itrack2 < pert->ntrack; itrack2++)
      for (int ixtrack2 = 0; ixtrack2 < pert->nxtrack; ixtrack2++)
	geo2cart(0, pert->lon[itrack2][ixtrack2],
		 pert->lat[itrack2][ixtrack2], xf[itrack2][ixtrack2]);

    /* Loop over storms... */
    for (size_t istorm = 0; istorm < nstorm; istorm++) {

      /* Loop along AIRS center track... */
      for (int itrack = 0; itrack < pert->ntrack; itrack++) {

	/* Get storm position... */
	if (get_storm_pos(nobs[istorm], time_wmo[istorm], lon_wmo[istorm],
			  lat_wmo[istorm], wind_wmo[istorm], pres_wmo[istorm],
			  pert->time[itrack][pert->nxtrack / 2], dt, st, xs,
			  &wind, &dwind, &pres, &dpres)) {

	  /* Get distance... */
	  double r2 = DIST2(xs, xf[itrack][pert->nxtrack / 2]);

	  /* Find best match... */
	  if (r2 < r2best) {

	    /* Save position... */
	    r2best = r2;
	    timebest = pert->time[itrack][pert->nxtrack / 2];
	    cart2geo(xs, &z, &lonbest, &latbest);

	    /* Save wind... */
	    windbest = wind;
	    dwindbest = dwind;
	    presbest = pres;
	    dpresbest = dpres;

	    /* Get BT data... */
	    n = 0;
	    bt8_min = 1e100;
	    bt4_mean = 0;
	    bt4_var = 0;
	    for (int itrack2 = GSL_MAX(itrack - ((int)(rmax / 17) + 1), 0);
		 itrack2 <= GSL_MIN(itrack + ((int) (rmax / 17) + 1),
				    pert->ntrack - 1); itrack2++)
	      for (int ixtrack2 = 0; ixtrack2 < pert->nxtrack; ixtrack2++) {

		/* Check data... */
		if (pert->time[itrack2][ixtrack2] < 0
		    || pert->lon[itrack2][ixtrack2] < -180
		    || pert->lon[itrack2][ixtrack2] > 180
		    || pert->lat[itrack2][ixtrack2] < -90
		    || pert->lat[itrack2][ixtrack2] > 90
		    || pert->pt[itrack2][ixtrack2] < -100
		    || pert->pt[itrack2][ixtrack2] > 100
		    || !gsl_finite(pert->bt[itrack2][ixtrack2])
		    || !gsl_finite(pert->pt[itrack2][ixtrack2])
		    || !gsl_finite(pert->var[itrack2][ixtrack2])
		    || !gsl_finite(pert->dc[itrack2][ixtrack2]))
		  continue;

		/* Check east/west filter... */
		lonsat = pert->lon[itrack2][ixtrack2];
		while (lonsat < 20)
		  lonsat += 360;
		lonstorm = lonbest;
		while (lonstorm < 20)
		  lonstorm += 360;
		if ((filter[0] == 'e' || filter[0] == 'E')
		    && lonsat < lonstorm)
		  continue;
		if ((filter[0] == 'w' || filter[0] == 'W')
		    && lonsat > lonstorm)
		  continue;

		/* Get distance... */
		if (DIST2(xs, xf[itrack2][ixtrack2]) < rmax * rmax) {
		  bt8_min = GSL_MIN(bt8_min, pert->dc[itrack2][ixtrack2]);
		  bt4_mean += pert->bt[itrack2][ixtrack2];
		  bt4_var += gsl_pow_2(pert->pt[itrack2][ixtrack2]);
		  n++;
		}
	      }
	  }
	}

	/* Output over poles... */
	if (fabs(pert->lat[itrack][pert->nxtrack / 2]) > 80.) {

	  /* Get and check ascending/descending flag... */
	  const int asc =
	    (pert->lat[itrack > 0 ? itrack : itrack + 1][pert->nxtrack / 2]
	     > pert->lat[itrack >
			 0 ? itrack - 1 : itrack][pert->nxtrack / 2]);
	  if ((set[0] == 'f' || set[0] == 'F')
	      || ((set[0] == 'a' || set[0] == 'A') && asc)
	      || ((set[0] == 'd' || set[0] == 'D') && !asc)) {

	    /* Check for match... */
	    if (r2best < 890. * 890.) {

	      /* Estimate noise... */
	      if (dt230 > 0) {
		const double nesr = NESR(230.0, dt230, nu);
		nedt = NEDT(bt4_mean / n, nesr, nu);
	      }
	      
	      /* Write output... */
	      if (n > 0)
		fprintf(out,
			"%lu %g %g %g %.2f %g %g %g %g %g %g %g %g %g %g %d\n",
			istorm, (timebest - time_wmo[istorm][0]) / 3600.,
			(timebest - time_max_wind[istorm]) / 3600.,
			(timebest - time_max_pres[istorm]) / 3600.,
			timebest, lonbest, latbest, sqrt(r2best), windbest,
			dwindbest, presbest, dpresbest, bt8_min, bt4_var / n,
			bt4_var / n - gsl_pow_2(nedt), n);
	    }
	  }

	  /* Reset... */
	  r2best = 1e100;
	}
      }
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

int get_storm_pos(
  int nobs,
  double time_wmo[NTIME],
  double lon_wmo[NTIME],
  double lat_wmo[NTIME],
  double wind_wmo[NTIME],
  double pres_wmo[NTIME],
  double t,
  int dt,
  int st,
  double x[3],
  double *wind,
  double *dwind,
  double *pres,
  double *dpres) {

  double x0[3], x1[3];

  /* Check time range... */
  if (t < time_wmo[0] || t > time_wmo[nobs - 1])
    return 0;

  /* Interpolate position... */
  int i = locate_irr(time_wmo, nobs, t);
  double w = (t - time_wmo[i]) / (time_wmo[i + 1] - time_wmo[i]);
  geo2cart(0, lon_wmo[i], lat_wmo[i], x0);
  geo2cart(0, lon_wmo[i + 1], lat_wmo[i + 1], x1);
  x[0] = (1 - w) * x0[0] + w * x1[0];
  x[1] = (1 - w) * x0[1] + w * x1[1];
  x[2] = (1 - w) * x0[2] + w * x1[2];

  /* Interpolate wind and pressure... */
  *pres = (1 - w) * pres_wmo[i] + w * pres_wmo[i + 1];
  *wind = (1 - w) * wind_wmo[i] + w * wind_wmo[i + 1];

  /* Get pressure and wind change... */
  *dpres = (pres_wmo[i + 1 + st] - pres_wmo[GSL_MAX(i - dt + st, 0)])
    / (time_wmo[i + 1 + st] - time_wmo[GSL_MAX(i - dt + st, 0)]) * 3600.;
  *dwind = (wind_wmo[i + 1 + st] - wind_wmo[GSL_MAX(i - dt + st, 0)])
    / (time_wmo[i + 1 + st] - time_wmo[GSL_MAX(i - dt + st, 0)]) * 3600.;

  return 1;
}

/*****************************************************************************/

void read_var(
  int ncid,
  const char varname[],
  size_t nstorm,
  int nobs[NSTORM],
  double x[NSTORM][NTIME]) {

  int varid;

  size_t count[2], start[2];

  /* Read pressure... */
  NC(nc_inq_varid(ncid, varname, &varid));
  for (size_t istorm = 0; istorm < nstorm; istorm++) {
    start[0] = istorm;
    start[1] = 0;
    count[0] = 1;
    count[1] = (size_t) nobs[istorm];
    NC(nc_get_vara_double(ncid, varid, start, count, x[istorm]));
  }
}
