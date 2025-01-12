#include "libairs.h"

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/* Number of 4 micron channels: */
#define N4 42

/* Number of 15 micron channels (low altitudes): */
#define N15_LOW 21

/* Number of 15 micron channels (high altitudes): */
#define N15_HIGH 2

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;

  static pert_t *pert_4mu, *pert_15mu_low, *pert_15mu_high;

  static wave_t wave;

  static double var_dh = 100.;

  static int list_4mu[N4]
    = { 2039, 2040, 2041, 2042, 2043, 2044, 2045, 2046, 2047, 2048,
    2049, 2050, 2051, 2052, 2053, 2054, 2055, 2056, 2057, 2058,
    2059, 2060, 2061, 2062, 2063, 2064, 2071, 2072, 2073, 2074,
    2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083, 2084,
    2085, 2086
  };

  static int list_15mu_low[N15_LOW]
    = { 4, 10, 16, 22, 29, 35, 41, 55, 83, 88, 94,
    100, 101, 106, 107, 112, 113, 118, 119, 124, 125
  };

  static int list_15mu_high[N15_HIGH]
  = { 74, 75 };

  static int ix, iy, dimid[2], i, n, ncid, track, track0, xtrack,
    time_varid, lon_varid, lat_varid, bt_4mu_varid, bt_4mu_pt_varid,
    bt_4mu_var_varid, bt_8mu_varid, bt_15mu_low_varid, bt_15mu_low_pt_varid,
    bt_15mu_low_var_varid, bt_15mu_high_varid, bt_15mu_high_pt_varid,
    bt_15mu_high_var_varid, iarg;

  static size_t start[2], count[2];

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <out.nc> <l1b_file1> [<l1b_file2> ...]");

  /* Allocate... */
  ALLOC(pert_4mu, pert_t, 1);
  ALLOC(pert_15mu_low, pert_t, 1);
  ALLOC(pert_15mu_high, pert_t, 1);

  /* ------------------------------------------------------------
     Read HDF files...
     ------------------------------------------------------------ */

  /* Loop over HDF files... */
  for (iarg = 2; iarg < argc; iarg++) {

    /* Read AIRS data... */
    printf("Read AIRS Level-1B data file: %s\n", argv[iarg]);
    airs_rad_rdr(argv[iarg], &airs_rad_gran);

    /* Flag bad observations... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
	for (i = 0; i < AIRS_RAD_CHANNEL; i++)
	  if ((airs_rad_gran.state[track][xtrack] != 0)
	      || (airs_rad_gran.ExcludedChans[i] > 2)
	      || (airs_rad_gran.CalChanSummary[i] & 8)
	      || (airs_rad_gran.CalChanSummary[i] & (32 + 64))
	      || (airs_rad_gran.CalFlag[track][i] & 16)
	      || (airs_rad_gran.Longitude[track][xtrack] < -180)
	      || (airs_rad_gran.Longitude[track][xtrack] > 180)
	      || (airs_rad_gran.Latitude[track][xtrack] < -90)
	      || (airs_rad_gran.Latitude[track][xtrack] > 90))
	    airs_rad_gran.radiances[track][xtrack][i] = GSL_NAN;
	  else
	    airs_rad_gran.radiances[track][xtrack][i] *= 0.001f;

    /* Save geolocation... */
    pert_4mu->ntrack += AIRS_RAD_GEOTRACK;
    if (pert_4mu->ntrack > PERT_NTRACK)
      ERRMSG("Too many granules!");
    pert_4mu->nxtrack = AIRS_RAD_GEOXTRACK;
    if (pert_4mu->nxtrack > PERT_NXTRACK)
      ERRMSG("Too many tracks!");
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
	pert_4mu->time[track0 + track][xtrack]
	  = airs_rad_gran.Time[track][xtrack] - 220838400.;
	pert_4mu->lon[track0 + track][xtrack]
	  = airs_rad_gran.Longitude[track][xtrack];
	pert_4mu->lat[track0 + track][xtrack]
	  = airs_rad_gran.Latitude[track][xtrack];
      }

    pert_15mu_low->ntrack += AIRS_RAD_GEOTRACK;
    if (pert_15mu_low->ntrack > PERT_NTRACK)
      ERRMSG("Too many granules!");
    pert_15mu_low->nxtrack = AIRS_RAD_GEOXTRACK;
    if (pert_15mu_low->nxtrack > PERT_NXTRACK)
      ERRMSG("Too many tracks!");
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
	pert_15mu_low->time[track0 + track][xtrack]
	  = airs_rad_gran.Time[track][xtrack] - 220838400.;
	pert_15mu_low->lon[track0 + track][xtrack]
	  = airs_rad_gran.Longitude[track][xtrack];
	pert_15mu_low->lat[track0 + track][xtrack]
	  = airs_rad_gran.Latitude[track][xtrack];
      }

    pert_15mu_high->ntrack += AIRS_RAD_GEOTRACK;
    if (pert_15mu_high->ntrack > PERT_NTRACK)
      ERRMSG("Too many granules!");
    pert_15mu_high->nxtrack = AIRS_RAD_GEOXTRACK;
    if (pert_15mu_high->nxtrack > PERT_NXTRACK)
      ERRMSG("Too many tracks!");
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
	pert_15mu_high->time[track0 + track][xtrack]
	  = airs_rad_gran.Time[track][xtrack] - 220838400.;
	pert_15mu_high->lon[track0 + track][xtrack]
	  = airs_rad_gran.Longitude[track][xtrack];
	pert_15mu_high->lat[track0 + track][xtrack]
	  = airs_rad_gran.Latitude[track][xtrack];
      }

    /* Get 8.1 micron brightness temperature... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
	pert_4mu->dc[track0 + track][xtrack]
	  = BRIGHT(airs_rad_gran.radiances[track][xtrack][1290],
		   airs_rad_gran.nominal_freq[1290]);

    /* Get 4.3 micron brightness temperature... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
	n = 0;
	for (i = 0; i < N4; i++)
	  if (gsl_finite(airs_rad_gran.radiances[track][xtrack][list_4mu[i]])) {
	    pert_4mu->bt[track0 + track][xtrack]
	      += BRIGHT(airs_rad_gran.radiances[track][xtrack][list_4mu[i]],
			airs_rad_gran.nominal_freq[list_4mu[i]]);
	    n++;
	  }
	if (n > 0.9 * N4)
	  pert_4mu->bt[track0 + track][xtrack] /= n;
	else
	  pert_4mu->bt[track0 + track][xtrack] = GSL_NAN;
      }

    /* Get 15 micron brightness temperature (low altitudes)... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
	n = 0;
	for (i = 0; i < N15_LOW; i++)
	  if (gsl_finite(airs_rad_gran.radiances
			 [track][xtrack][list_15mu_low[i]])) {
	    pert_15mu_low->bt[track0 + track][xtrack]
	      +=
	      BRIGHT(airs_rad_gran.radiances[track][xtrack][list_15mu_low[i]],
		     airs_rad_gran.nominal_freq[list_15mu_low[i]]);
	    n++;
	  }
	if (n > 0.9 * N15_LOW)
	  pert_15mu_low->bt[track0 + track][xtrack] /= n;
	else
	  pert_15mu_low->bt[track0 + track][xtrack] = GSL_NAN;
      }

    /* Get 15 micron brightness temperature (high altitudes)... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
	n = 0;
	for (i = 0; i < N15_HIGH; i++)
	  if (gsl_finite(airs_rad_gran.radiances
			 [track][xtrack][list_15mu_high[i]])) {
	    pert_15mu_high->bt[track0 + track][xtrack]
	      += BRIGHT(airs_rad_gran.radiances
			[track][xtrack][list_15mu_high[i]],
			airs_rad_gran.nominal_freq[list_15mu_high[i]]);
	    n++;
	  }
	if (n > 0.9 * N15_HIGH)
	  pert_15mu_high->bt[track0 + track][xtrack] /= n;
	else
	  pert_15mu_high->bt[track0 + track][xtrack] = GSL_NAN;
      }

    /* Increment track counter... */
    track0 += AIRS_RAD_GEOTRACK;
  }

  /* ------------------------------------------------------------
     Calculate perturbations and variances...
     ------------------------------------------------------------ */

  /* Convert to wave analysis struct... */
  pert2wave(pert_4mu, &wave,
	    0, pert_4mu->ntrack - 1, 0, pert_4mu->nxtrack - 1);

  /* Estimate background... */
  background_poly(&wave, 5, 0);

  /* Compute variance... */
  variance(&wave, var_dh);

  /* Copy data... */
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++) {
      pert_4mu->pt[iy][ix] = wave.pt[ix][iy];
      pert_4mu->var[iy][ix] = wave.var[ix][iy];
    }

  /* Convert to wave analysis struct... */
  pert2wave(pert_15mu_low, &wave,
	    0, pert_15mu_low->ntrack - 1, 0, pert_15mu_low->nxtrack - 1);

  /* Estimate background... */
  background_poly(&wave, 5, 0);

  /* Compute variance... */
  variance(&wave, var_dh);

  /* Copy data... */
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++) {
      pert_15mu_low->pt[iy][ix] = wave.pt[ix][iy];
      pert_15mu_low->var[iy][ix] = wave.var[ix][iy];
    }

  /* Convert to wave analysis struct... */
  pert2wave(pert_15mu_high, &wave,
	    0, pert_15mu_high->ntrack - 1, 0, pert_15mu_high->nxtrack - 1);

  /* Estimate background... */
  background_poly(&wave, 5, 0);

  /* Compute variance... */
  variance(&wave, var_dh);

  /* Copy data... */
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++) {
      pert_15mu_high->pt[iy][ix] = wave.pt[ix][iy];
      pert_15mu_high->var[iy][ix] = wave.var[ix][iy];
    }

  /* ------------------------------------------------------------
     Write to netCDF file...
     ------------------------------------------------------------ */

  /* Create netCDF file... */
  NC(nc_create(argv[1], NC_CLOBBER, &ncid));

  /* Set dimensions... */
  NC(nc_def_dim(ncid, "NTRACK", NC_UNLIMITED, &dimid[0]));
  NC(nc_def_dim(ncid, "NXTRACK", AIRS_RAD_GEOXTRACK, &dimid[1]));

  /* Add variables... */
  NC(nc_def_var(ncid, "time", NC_DOUBLE, 2, dimid, &time_varid));
  add_att(ncid, time_varid, "s", "time (seconds since 2000-01-01T00:00Z)");
  NC(nc_def_var(ncid, "lon", NC_DOUBLE, 2, dimid, &lon_varid));
  add_att(ncid, lon_varid, "deg", "footprint longitude");
  NC(nc_def_var(ncid, "lat", NC_DOUBLE, 2, dimid, &lat_varid));
  add_att(ncid, lat_varid, "deg", "footprint latitude");

  NC(nc_def_var(ncid, "bt_8mu", NC_FLOAT, 2, dimid, &bt_8mu_varid));
  add_att(ncid, bt_8mu_varid, "K", "brightness temperature at 8.1 micron");

  NC(nc_def_var(ncid, "bt_4mu", NC_FLOAT, 2, dimid, &bt_4mu_varid));
  add_att(ncid, bt_4mu_varid, "K", "brightness temperature" " at 4.3 micron");
  NC(nc_def_var(ncid, "bt_4mu_pt", NC_FLOAT, 2, dimid, &bt_4mu_pt_varid));
  add_att(ncid, bt_4mu_pt_varid, "K", "brightness temperature perturbation"
	  " at 4.3 micron");
  NC(nc_def_var(ncid, "bt_4mu_var", NC_FLOAT, 2, dimid, &bt_4mu_var_varid));
  add_att(ncid, bt_4mu_var_varid, "K^2", "brightness temperature variance"
	  " at 4.3 micron");

  NC(nc_def_var(ncid, "bt_15mu_low", NC_FLOAT, 2, dimid, &bt_15mu_low_varid));
  add_att(ncid, bt_15mu_low_varid, "K", "brightness temperature"
	  " at 15 micron (low altitudes)");
  NC(nc_def_var(ncid, "bt_15mu_low_pt", NC_FLOAT, 2, dimid,
		&bt_15mu_low_pt_varid));
  add_att(ncid, bt_15mu_low_pt_varid, "K",
	  "brightness temperature perturbation"
	  " at 15 micron (low altitudes)");
  NC(nc_def_var
     (ncid, "bt_15mu_low_var", NC_FLOAT, 2, dimid, &bt_15mu_low_var_varid));
  add_att(ncid, bt_15mu_low_var_varid, "K^2",
	  "brightness temperature variance" " at 15 micron (low altitudes)");

  NC(nc_def_var(ncid, "bt_15mu_high", NC_FLOAT, 2, dimid,
		&bt_15mu_high_varid));
  add_att(ncid, bt_15mu_high_varid, "K", "brightness temperature"
	  " at 15 micron (high altitudes)");
  NC(nc_def_var(ncid, "bt_15mu_high_pt", NC_FLOAT, 2, dimid,
		&bt_15mu_high_pt_varid));
  add_att(ncid, bt_15mu_high_pt_varid, "K",
	  "brightness temperature perturbation"
	  " at 15 micron (high altitudes)");
  NC(nc_def_var
     (ncid, "bt_15mu_high_var", NC_FLOAT, 2, dimid, &bt_15mu_high_var_varid));
  add_att(ncid, bt_15mu_high_var_varid, "K^2",
	  "brightness temperature variance" " at 15 micron (high altitudes)");

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Loop over tracks... */
  for (track = 0; track < pert_4mu->ntrack; track++) {

    /* Set array sizes... */
    start[0] = (size_t) track;
    start[1] = 0;
    count[0] = 1;
    count[1] = (size_t) pert_4mu->nxtrack;

    /* Write data... */
    NC(nc_put_vara_double(ncid, time_varid, start, count,
			  pert_4mu->time[track]));
    NC(nc_put_vara_double(ncid, lon_varid, start, count,
			  pert_4mu->lon[track]));
    NC(nc_put_vara_double(ncid, lat_varid, start, count,
			  pert_4mu->lat[track]));

    NC(nc_put_vara_double(ncid, bt_8mu_varid, start, count,
			  pert_4mu->dc[track]));

    NC(nc_put_vara_double(ncid, bt_4mu_varid, start, count,
			  pert_4mu->bt[track]));
    NC(nc_put_vara_double(ncid, bt_4mu_pt_varid, start, count,
			  pert_4mu->pt[track]));
    NC(nc_put_vara_double(ncid, bt_4mu_var_varid, start, count,
			  pert_4mu->var[track]));

    NC(nc_put_vara_double(ncid, bt_15mu_low_varid, start, count,
			  pert_15mu_low->bt[track]));
    NC(nc_put_vara_double(ncid, bt_15mu_low_pt_varid, start, count,
			  pert_15mu_low->pt[track]));
    NC(nc_put_vara_double(ncid, bt_15mu_low_var_varid, start, count,
			  pert_15mu_low->var[track]));

    NC(nc_put_vara_double(ncid, bt_15mu_high_varid, start, count,
			  pert_15mu_high->bt[track]));
    NC(nc_put_vara_double(ncid, bt_15mu_high_pt_varid, start, count,
			  pert_15mu_high->pt[track]));
    NC(nc_put_vara_double(ncid, bt_15mu_high_var_varid, start, count,
			  pert_15mu_high->var[track]));
  }

  /* Close file... */
  NC(nc_close(ncid));

  /* Free... */
  free(pert_4mu);
  free(pert_15mu_low);
  free(pert_15mu_high);

  return EXIT_SUCCESS;
}
