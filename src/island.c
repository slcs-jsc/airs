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
  Analyse gravity wave data for remote islands.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  static wave_t *wave;

  static FILE *in, *out;

  static char pertname[LEN], ncfile[LEN];

  static double gauss_fwhm, var_dh, orblat, lon0, lat0, dlon, dlat, offset,
    ebt, emu, enoise, evar, wbt, wmu, wnoise, wvar, etime, wtime,
    dt230, nu, nesr, aux;

  static int iarg, ix, iy, itrack, itrack2, ixtrack, bg_poly_x, bg_poly_y,
    bg_smooth_x, bg_smooth_y, orb, orb_old = -1, en, wn, ncid, dimid[2],
    time_varid, track_varid, np_east_varid, var_east_varid,
    np_west_varid, var_west_varid, year_varid, doy_varid,
    track, year, mon, day, doy, iaux;

  static size_t count[2] = { 1, 1 }, start[2];

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <var.tab> <pert1.nc> [<pert2.nc> ...]");

  /* Get control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  lon0 = scan_ctl(argc, argv, "LON0", -1, "", NULL);
  lat0 = scan_ctl(argc, argv, "LAT0", -1, "", NULL);
  dlon = scan_ctl(argc, argv, "DLON", -1, "", NULL);
  dlat = scan_ctl(argc, argv, "DLAT", -1, "", NULL);
  offset = scan_ctl(argc, argv, "OFFSET", -1, "1", NULL);
  bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "0", NULL);
  bg_poly_y = (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  bg_smooth_x = (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  bg_smooth_y = (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "0", NULL);
  gauss_fwhm = scan_ctl(argc, argv, "GAUSS_FWHM", -1, "0", NULL);
  var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "0", NULL);
  orblat = scan_ctl(argc, argv, "ORBLAT", -1, "0", NULL);
  dt230 = scan_ctl(argc, argv, "DT230", -1, "0.16", NULL);
  nu = scan_ctl(argc, argv, "NU", -1, "2345.0", NULL);
  scan_ctl(argc, argv, "NCFILE", -1, "-", ncfile);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Create file... */
  printf("Write variance statistics: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time [s]\n"
	  "# $2  = orbit number\n"
	  "# $3  = eastern box: number of footprints\n"
	  "# $4  = eastern box: variance [K^2]\n"
	  "# $5  = eastern box: mean background temperature [K]\n"
	  "# $6  = eastern box: noise estimate [K]\n"
	  "# $7  = western box: number of footprints\n"
	  "# $8  = western box: variance [K^2]\n"
	  "# $9  = western box: mean background temperature [K]\n"
	  "# $10 = western box: noise estimate [K]\n\n");

  /* Create netCDF file... */
  if (ncfile[0] != '-') {

    /* Create file... */
    printf("Write variance statistics: %s\n", ncfile);
    NC(nc_create(ncfile, NC_CLOBBER, &ncid));

    /* Set dimensions... */
    NC(nc_def_dim(ncid, "NP", NC_UNLIMITED, &dimid[0]));

    /* Add attributes... */
    aux = lon0;
    nc_put_att_double(ncid, NC_GLOBAL, "box_east_lon0", NC_DOUBLE, 1, &aux);
    aux = lon0 + dlon;
    nc_put_att_double(ncid, NC_GLOBAL, "box_east_lon1", NC_DOUBLE, 1, &aux);
    aux = lat0 - 0.5 * dlat;
    nc_put_att_double(ncid, NC_GLOBAL, "box_east_lat0", NC_DOUBLE, 1, &aux);
    aux = lat0 + 0.5 * dlat;
    nc_put_att_double(ncid, NC_GLOBAL, "box_east_lat1", NC_DOUBLE, 1, &aux);
    aux = lon0 - dlon - offset;
    nc_put_att_double(ncid, NC_GLOBAL, "box_west_lon0", NC_DOUBLE, 1, &aux);
    aux = lon0 - offset;
    nc_put_att_double(ncid, NC_GLOBAL, "box_west_lon1", NC_DOUBLE, 1, &aux);
    aux = lat0 - 0.5 * dlat;
    nc_put_att_double(ncid, NC_GLOBAL, "box_west_lat0", NC_DOUBLE, 1, &aux);
    aux = lat0 + 0.5 * dlat;
    nc_put_att_double(ncid, NC_GLOBAL, "box_west_lat1", NC_DOUBLE, 1, &aux);

    /* Add variables... */
    NC(nc_def_var(ncid, "time", NC_DOUBLE, 1, dimid, &time_varid));
    add_att(ncid, time_varid, "s", "time (seconds since 2000-01-01T00:00Z)");
    NC(nc_def_var(ncid, "year", NC_INT, 1, dimid, &year_varid));
    add_att(ncid, year_varid, "1", "year");
    NC(nc_def_var(ncid, "doy", NC_INT, 1, dimid, &doy_varid));
    add_att(ncid, doy_varid, "1", "day of year");
    NC(nc_def_var(ncid, "track", NC_INT, 1, dimid, &track_varid));
    add_att(ncid, track_varid, "1", "along-track index");
    NC(nc_def_var(ncid, "var_east", NC_DOUBLE, 1, dimid, &var_east_varid));
    add_att(ncid, var_east_varid, "K^2", "BT variance (east)");
    NC(nc_def_var(ncid, "var_west", NC_DOUBLE, 1, dimid, &var_west_varid));
    add_att(ncid, var_west_varid, "K^2", "BT variance (west)");
    NC(nc_def_var(ncid, "np_east", NC_INT, 1, dimid, &np_east_varid));
    add_att(ncid, np_east_varid, "1", "number of footprints (east)");
    NC(nc_def_var(ncid, "np_west", NC_INT, 1, dimid, &np_west_varid));
    add_att(ncid, np_west_varid, "1", "number of footprints (west)");

    /* Leave define mode... */
    NC(nc_enddef(ncid));
  }

  /* Loop over perturbation files... */
  for (iarg = 3; iarg < argc; iarg++) {

    /* Check filename... */
    if (!strcmp(argv[iarg], ncfile))
      continue;

    /* Initialize... */
    orb = 0;

    /* Read perturbation data... */
    if (!(in = fopen(argv[iarg], "r")))
      continue;
    else {
      fclose(in);
      read_pert(argv[iarg], pertname, pert);
    }

    /* Recalculate background and perturbations... */
    if (bg_poly_x > 0 || bg_poly_y > 0 ||
	bg_smooth_x > 0 || bg_smooth_y > 0 || gauss_fwhm > 0 || var_dh > 0) {

      /* Allocate... */
      ALLOC(wave, wave_t, 1);

      /* Convert to wave analysis struct... */
      pert2wave(pert, wave, 0, pert->ntrack - 1, 0, pert->nxtrack - 1);

      /* Estimate background... */
      background_poly(wave, bg_poly_x, bg_poly_y);
      background_smooth(wave, bg_smooth_x, bg_smooth_y);

      /* Gaussian filter... */
      gauss(wave, gauss_fwhm);

      /* Compute variance... */
      variance(wave, var_dh);

      /* Copy data... */
      for (ix = 0; ix < wave->nx; ix++)
	for (iy = 0; iy < wave->ny; iy++) {
	  pert->pt[iy][ix] = wave->pt[ix][iy];
	  pert->var[iy][ix] = wave->var[ix][iy];
	}

      /* Free... */
      free(wave);
    }

    /* Detection... */
    for (itrack = 0; itrack < pert->ntrack; itrack++)
      for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {

	/* Check data... */
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

	/* Count orbits... */
	if (itrack > 0 && ixtrack == pert->nxtrack / 2)
	  if (pert->lat[itrack - 1][ixtrack] <= orblat
	      && pert->lat[itrack][ixtrack] >= orblat)
	    orb++;
	if (orb != orb_old) {

	  /* Set orbit index... */
	  orb_old = orb;

	  /* Write output... */
	  if (en > 0 && wn > 0) {

	    /* Estimate noise... */
	    if (dt230 > 0) {
	      nesr = PLANCK(230.0 + dt230, nu) - PLANCK(230.0, nu);
	      enoise = BRIGHT(PLANCK(ebt / en, nu) + nesr, nu) - ebt / en;
	      wnoise = BRIGHT(PLANCK(wbt / wn, nu) + nesr, nu) - wbt / wn;
	    }

	    /* Write output... */
	    fprintf(out, "%.2f %d %d %g %g %g %d %g %g %g\n", etime / en, orb,
		    en, evar / en - gsl_pow_2(emu / en), ebt / en, enoise,
		    wn, wvar / wn - gsl_pow_2(wmu / wn), wbt / wn, wnoise);

	    /* Write to netCDF file... */
	    if (ncfile[0] != '-') {

	      /* Get year and doy... */
	      jsec2time(etime / en, &year, &mon, &day, &iaux, &iaux, &iaux,
			&aux);
	      day2doy(year, mon, day, &doy);

	      /* Find along-track index... */
	      track = 0;
	      for (itrack2 = 0; itrack2 < pert->ntrack; itrack2++)
		if (fabs(pert->time[itrack2][0] - etime / en)
		    < fabs(pert->time[track][0] - etime / en))
		  track = itrack2;

	      /* Write data... */
	      aux = etime / en;
	      NC(nc_put_vara_double(ncid, time_varid, start, count, &aux));
	      NC(nc_put_vara_int(ncid, year_varid, start, count, &year));
	      NC(nc_put_vara_int(ncid, doy_varid, start, count, &doy));
	      NC(nc_put_vara_int(ncid, track_varid, start, count, &track));
	      NC(nc_put_vara_int(ncid, np_east_varid, start, count, &en));
	      aux = evar / en - gsl_pow_2(emu / en) - gsl_pow_2(enoise);
	      NC(nc_put_vara_double
		 (ncid, var_east_varid, start, count, &aux));
	      NC(nc_put_vara_int(ncid, np_west_varid, start, count, &wn));
	      aux = wvar / wn - gsl_pow_2(wmu / wn) - gsl_pow_2(wnoise);
	      NC(nc_put_vara_double
		 (ncid, var_west_varid, start, count, &aux));

	      /* Increment data point counter... */
	      start[0]++;
	    }
	  }

	  /* Initialize... */
	  etime = wtime = 0;
	  evar = wvar = 0;
	  emu = wmu = 0;
	  ebt = wbt = 0;
	  en = wn = 0;
	}

	/* Check if footprint is in eastern box... */
	if (pert->lon[itrack][ixtrack] >= lon0
	    && pert->lon[itrack][ixtrack] <= lon0 + dlon
	    && pert->lat[itrack][ixtrack] >= lat0 - dlat / 2.
	    && pert->lat[itrack][ixtrack] <= lat0 + dlat / 2.) {

	  etime += pert->time[itrack][ixtrack];
	  emu += pert->pt[itrack][ixtrack];
	  evar += gsl_pow_2(pert->pt[itrack][ixtrack]);
	  ebt += pert->bt[itrack][ixtrack];
	  en++;
	}

	/* Check if footprint is in western box... */
	if (pert->lon[itrack][ixtrack] >= lon0 - offset - dlon
	    && pert->lon[itrack][ixtrack] <= lon0 - offset
	    && pert->lat[itrack][ixtrack] >= lat0 - dlat / 2.
	    && pert->lat[itrack][ixtrack] <= lat0 + dlat / 2.) {

	  wtime += pert->time[itrack][ixtrack];
	  wmu += pert->pt[itrack][ixtrack];
	  wvar += gsl_pow_2(pert->pt[itrack][ixtrack]);
	  wbt += pert->bt[itrack][ixtrack];
	  wn++;
	}
      }

    /* Write output for last orbit... */
    if (en > 0 && wn > 0) {

      /* Estimate noise... */
      if (dt230 > 0) {
	nesr = PLANCK(230.0 + dt230, nu) - PLANCK(230.0, nu);
	enoise = BRIGHT(PLANCK(ebt / en, nu) + nesr, nu) - ebt / en;
	wnoise = BRIGHT(PLANCK(wbt / wn, nu) + nesr, nu) - wbt / wn;
      }

      /* Write output... */
      fprintf(out, "%.2f %d %d %g %g %g %d %g %g %g\n", etime / en, orb,
	      en, evar / en - gsl_pow_2(emu / en), ebt / en, enoise,
	      wn, wvar / wn - gsl_pow_2(wmu / wn), wbt / wn, wnoise);

      /* Write to netCDF file... */
      if (ncfile[0] != '-') {

	/* Get year and doy... */
	jsec2time(etime / en, &year, &mon, &day, &iaux, &iaux, &iaux, &aux);
	day2doy(year, mon, day, &doy);

	/* Find along-track index... */
	track = 0;
	for (itrack2 = 0; itrack2 < pert->ntrack; itrack2++)
	  if (fabs(pert->time[itrack2][0] - etime / en)
	      < fabs(pert->time[track][0] - etime / en))
	    track = itrack2;

	/* Write data... */
	aux = etime / en;
	NC(nc_put_vara_double(ncid, time_varid, start, count, &aux));
	NC(nc_put_vara_int(ncid, year_varid, start, count, &year));
	NC(nc_put_vara_int(ncid, doy_varid, start, count, &doy));
	NC(nc_put_vara_int(ncid, track_varid, start, count, &track));
	NC(nc_put_vara_int(ncid, np_east_varid, start, count, &en));
	aux = evar / en - gsl_pow_2(emu / en) - gsl_pow_2(enoise);
	NC(nc_put_vara_double(ncid, var_east_varid, start, count, &aux));
	NC(nc_put_vara_int(ncid, np_west_varid, start, count, &wn));
	aux = wvar / wn - gsl_pow_2(wmu / wn) - gsl_pow_2(wnoise);
	NC(nc_put_vara_double(ncid, var_west_varid, start, count, &aux));

	/* Increment data point counter... */
	start[0]++;
      }
    }
  }

  /* Close file... */
  fclose(out);

  /* Close file... */
  if (ncfile[0] != '-')
    NC(nc_close(ncid));

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}
