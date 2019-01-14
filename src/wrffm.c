#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum WRF dimensions. */
#define NLON 450
#define NLAT 450
#define NZ 150

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Estimate background... */
void background(
  double temp[NLAT][NLON],
  double pt[NLAT][NLON],
  int nlat,
  int nlon,
  int dlat,
  int dlon);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static atm_t atm;
  static obs_t obs;

  static double z[NZ][NLAT][NLON], z0[NZ][NLAT][NLON],
    p[NZ][NLAT][NLON], p0[NZ][NLAT][NLON], t[NZ][NLAT][NLON],
    lon[NLAT][NLON], lat[NLAT][NLON], temp[NLAT][NLON], pt[NLAT][NLON],
    x0[3], x1[NLAT][NLON][3], w, wsum, rmax2 = 50. * 50., fwhm = 20., var_dh;

  static int id, ix, iy, oit, ncid, dimid, varid, ilon, ilat, iz,
    ncrop, nlon, nlat, nz, nz2, ntime;

  static size_t start[10], count[10], rs;

  wave_t *wave_airs, *wave_wrf;

  FILE *out;

  /* ------------------------------------------------------------
     Get control parameters...
     ------------------------------------------------------------ */

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <wrf.nc> <it> <wrf.tab> <rad.tab> "
	   "<wave_airs.tab> <wave_wrf.tab>");

  /* Get arguments... */
  oit = atoi(argv[3]);

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Set control parameters... */
  ctl.write_bbt = 1;

  /* Get control parameters... */
  ncrop = (int) scan_ctl(argc, argv, "NCROP", -1, "0", NULL);
  var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "0", NULL);

  /* Allocate... */
  ALLOC(wave_airs, wave_t, 1);
  ALLOC(wave_wrf, wave_t, 1);

  /* ------------------------------------------------------------
     Read WRF data...
     ------------------------------------------------------------ */

  /* Open file... */
  printf("Read WRF data: %s\n", argv[2]);
  NC(nc_open(argv[2], NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "Time", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  ntime = (int) rs;
  if (oit >= ntime)
    ERRMSG("Timestep out of range!");

  NC(nc_inq_dimid(ncid, "bottom_top", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nz = (int) rs;
  if (nz > NZ)
    ERRMSG("Too many altitudes!");

  NC(nc_inq_dimid(ncid, "bottom_top_stag", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nz2 = (int) rs;
  if (nz2 > NZ)
    ERRMSG("Too many altitudes!");

  NC(nc_inq_dimid(ncid, "south_north", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nlat = (int) rs;
  if (nlat > NLAT)
    ERRMSG("Too many latitudes!");

  NC(nc_inq_dimid(ncid, "west_east", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nlon = (int) rs;
  if (nlon > NLON)
    ERRMSG("Too many longitudes!");

  /* Read latitudes... */
  NC(nc_inq_varid(ncid, "XLAT", &varid));
  for (ilat = 0; ilat < nlat; ilat++) {
    start[0] = (size_t) oit;
    start[1] = (size_t) ilat;
    start[2] = 0;
    count[0] = 1;
    count[1] = 1;
    count[2] = (size_t) nlon;
    NC(nc_get_vara_double(ncid, varid, start, count, lat[ilat]));
  }

  /* Read longitudes... */
  NC(nc_inq_varid(ncid, "XLONG", &varid));
  for (ilat = 0; ilat < nlat; ilat++) {
    start[0] = (size_t) oit;
    start[1] = (size_t) ilat;
    start[2] = 0;
    count[0] = 1;
    count[1] = 1;
    count[2] = (size_t) nlon;
    NC(nc_get_vara_double(ncid, varid, start, count, lon[ilat]));
  }

  /* Read theta perturbation... */
  NC(nc_inq_varid(ncid, "T", &varid));
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++) {
      start[0] = (size_t) oit;
      start[1] = (size_t) iz;
      start[2] = (size_t) ilat;
      start[3] = 0;
      count[0] = 1;
      count[1] = 1;
      count[2] = 1;
      count[3] = (size_t) nlon;
      NC(nc_get_vara_double(ncid, varid, start, count, t[iz][ilat]));
    }

  /* Read geopotential perturbation... */
  NC(nc_inq_varid(ncid, "PH", &varid));
  for (iz = 0; iz < nz2; iz++)
    for (ilat = 0; ilat < nlat; ilat++) {
      start[0] = (size_t) oit;
      start[1] = (size_t) iz;
      start[2] = (size_t) ilat;
      start[3] = 0;
      count[0] = 1;
      count[1] = 1;
      count[2] = 1;
      count[3] = (size_t) nlon;
      NC(nc_get_vara_double(ncid, varid, start, count, z[iz][ilat]));
    }

  /* Read geopotential base... */
  NC(nc_inq_varid(ncid, "PHB", &varid));
  for (iz = 0; iz < nz2; iz++)
    for (ilat = 0; ilat < nlat; ilat++) {
      start[0] = (size_t) oit;
      start[1] = (size_t) iz;
      start[2] = (size_t) ilat;
      start[3] = 0;
      count[0] = 1;
      count[1] = 1;
      count[2] = 1;
      count[3] = (size_t) nlon;
      NC(nc_get_vara_double(ncid, varid, start, count, z0[iz][ilat]));
    }

  /* Read pressure perturbation... */
  NC(nc_inq_varid(ncid, "P", &varid));
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++) {
      start[0] = (size_t) oit;
      start[1] = (size_t) iz;
      start[2] = (size_t) ilat;
      start[3] = 0;
      count[0] = 1;
      count[1] = 1;
      count[2] = 1;
      count[3] = (size_t) nlon;
      NC(nc_get_vara_double(ncid, varid, start, count, p[iz][ilat]));
    }

  /* Read pressure base... */
  NC(nc_inq_varid(ncid, "PB", &varid));
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++) {
      start[0] = (size_t) oit;
      start[1] = (size_t) iz;
      start[2] = (size_t) ilat;
      start[3] = 0;
      count[0] = 1;
      count[1] = 1;
      count[2] = 1;
      count[3] = (size_t) nlon;
      NC(nc_get_vara_double(ncid, varid, start, count, p0[iz][ilat]));
    }

  /* Close file... */
  NC(nc_close(ncid));

  /* ------------------------------------------------------------
     Convert WRF data...
     ------------------------------------------------------------ */

  /* Adjust longitudes... */
  for (ilat = 0; ilat < nlat; ilat++)
    for (ilon = 0; ilon < nlon; ilon++)
      if (lon[ilat][ilon] > 180)
	lon[ilat][ilon] -= 360;

  /* Get altitudes... */
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	z[iz][ilat][ilon]
	  = 0.5 * (z[iz + 1][ilat][ilon] + z0[iz + 1][ilat][ilon]
		   + z[iz][ilat][ilon] + z0[iz][ilat][ilon]) / G0 / 1000.;

  /* Get pressure... */
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	p[iz][ilat][ilon]
	  = (p[iz][ilat][ilon] + p0[iz][ilat][ilon]) / 100.;

  /* Get temperature... */
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	t[iz][ilat][ilon]
	  = (t[iz][ilat][ilon] + 300.) / pow(1000. / p[iz][ilat][ilon],
					     0.286);

  /* ------------------------------------------------------------
     Write WRF data to ASCII...
     ------------------------------------------------------------ */

  /* Check filename... */
  if (argv[4][0] != '-') {

    /* Create file... */
    printf("Write WRF data: %s\n", argv[4]);
    if (!(out = fopen(argv[4], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = altitude index\n"
	    "# $2 = altitude [km]\n"
	    "# $3 = longitude [deg]\n"
	    "# $4 = latitude [deg]\n"
	    "# $5 = pressure [hPa]\n" "# $6 = temperature [K]\n");

    /* Write output... */
    for (iz = 0; iz < nz; iz++)
      for (ilon = 0; ilon < nlon; ilon++) {
	fprintf(out, "\n");
	for (ilat = 0; ilat < nlat; ilat++)
	  fprintf(out, "%d %g %g %g %g %g\n", iz, z[iz][ilat][ilon],
		  lon[ilat][ilon], lat[ilat][ilon],
		  p[iz][ilat][ilon], t[iz][ilat][ilon]);
      }

    /* Close file... */
    fclose(out);
  }

  /* ------------------------------------------------------------
     Run forward model...
     ------------------------------------------------------------ */

  /* Loop over latitudes... */
  for (ilat = 0; ilat < nlat; ilat++) {

    /* Write info... */
    printf("  Compute latitude %d / %d ...\n", ilat + 1, nlat);

    /* Loop over longitudes... */
    for (ilon = 0; ilon < nlon; ilon++) {

      /* Set altitude levels... */
      atm.np = 0;
      for (iz = 0; iz < nz; iz++)
	if (gsl_finite(gsl_finite(t[iz][ilat][ilon]))
	    && t[iz][ilat][ilon] > 100 && t[iz][ilat][ilon] < 400
	    && z[iz][ilat][ilon] > 10) {
	  atm.z[atm.np] = z[iz][ilat][ilon];
	  if ((++atm.np) >= NP)
	    ERRMSG("Too many altitudes!");
	}

      /* Add top level... */
      atm.z[atm.np] = 90.;
      if ((++atm.np) >= NP)
	ERRMSG("Too many altitudes!");

      /* Initialize with climatological data... */
      climatology(&ctl, &atm);

      /* Set temperature and pressure... */
      atm.np = 0;
      for (iz = 0; iz < nz; iz++)
	if (gsl_finite(t[iz][ilat][ilon])
	    && t[iz][ilat][ilon] > 100 && t[iz][ilat][ilon] < 400
	    && z[iz][ilat][ilon] > 10) {
	  atm.p[atm.np] = p[iz][ilat][ilon];
	  atm.t[atm.np] = t[iz][ilat][ilon];
	  atm.np++;
	}

      /* Add top level... */
      atm.np++;

      /* Set observation data... */
      obs.nr = 1;
      obs.obsz[0] = 700;

      /* Run forward model... */
      formod(&ctl, &atm, &obs);

      /* Get mean brightness temperature... */
      temp[ilat][ilon] = 0;
      for (id = 0; id < ctl.nd; id++)
	temp[ilat][ilon] += obs.rad[id][0] / ctl.nd;
    }
  }

  /* Crop at boundaries... */
  for (ilat = 0; ilat < ncrop; ilat++)
    for (ilon = 0; ilon < nlon; ilon++)
      temp[ilat][ilon] = GSL_NAN;
  for (ilat = nlat - ncrop; ilat < nlat; ilat++)
    for (ilon = 0; ilon < nlon; ilon++)
      temp[ilat][ilon] = GSL_NAN;
  for (ilon = 0; ilon < ncrop; ilon++)
    for (ilat = 0; ilat < nlat; ilat++)
      temp[ilat][ilon] = GSL_NAN;
  for (ilon = nlon - ncrop; ilon < nlon; ilon++)
    for (ilat = 0; ilat < nlat; ilat++)
      temp[ilat][ilon] = GSL_NAN;

  /* Get perturbations... */
  background(temp, pt, nlat, nlon, 10, 10);

  /* ------------------------------------------------------------
     Save forward model output...
     ------------------------------------------------------------ */

  /* Check filename... */
  if (argv[5][0] != '-') {

    /* Create file... */
    printf("Write radiance data: %s\n", argv[5]);
    if (!(out = fopen(argv[5], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = longitude [deg]\n"
	    "# $2 = latitude [deg]\n"
	    "# $3 = WRF brightness temperature [K]\n"
	    "# $4 = WRF brightness temperature perturbation [K]\n");

    /* Write output... */
    for (ilat = 0; ilat < nlat; ilat++) {
      fprintf(out, "\n");
      for (ilon = 0; ilon < nlon; ilon++)
	fprintf(out, "%g %g %g %g\n", lon[ilat][ilon], lat[ilat][ilon],
		temp[ilat][ilon], pt[ilat][ilon]);
    }

    /* Close file... */
    fclose(out);
  }

  /* ------------------------------------------------------------
     Read AIRS radiance map and resample model data...
     ------------------------------------------------------------ */

  /* Check filename... */
  if (argv[6][0] != '-') {

    /* Read AIRS wave file... */
    read_wave(argv[6], wave_airs);
    memcpy(wave_wrf, wave_airs, sizeof(wave_t));

    /* Get Cartesian coordinates for model grid... */
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	geo2cart(0, lon[ilat][ilon], lat[ilat][ilon], x1[ilat][ilon]);

    /* Loop over AIRS geolocations... */
    for (ix = 0; ix < wave_airs->nx; ix++)
      for (iy = 0; iy < wave_airs->ny; iy++) {

	/* Write info... */
	if (iy == 0)
	  printf("  Average for xtrack %d / %d ...\n", ix + 1, wave_airs->nx);

	/* Init... */
	wsum = 0;
	wave_wrf->temp[ix][iy] = 0;
	wave_wrf->bg[ix][iy] = 0;
	wave_wrf->pt[ix][iy] = 0;
	wave_wrf->var[ix][iy] = 0;

	/* Average... */
	geo2cart(0, wave_airs->lon[ix][iy], wave_airs->lat[ix][iy], x0);
	for (ilat = 0; ilat < nlat; ilat++)
	  for (ilon = 0; ilon < nlon; ilon++)
	    if (DIST2(x0, x1[ilat][ilon]) <= rmax2) {
	      w =
		exp(-DIST2(x0, x1[ilat][ilon]) /
		    (2. * gsl_pow_2(fwhm / 2.3548)));
	      wave_wrf->temp[ix][iy] += w * temp[ilat][ilon];
	      wave_wrf->bg[ix][iy] += w * (temp[ilat][ilon] - pt[ilat][ilon]);
	      wave_wrf->pt[ix][iy] += w * pt[ilat][ilon];
	      wsum += w;
	    }

	/* Normalize... */
	if (wsum > 0) {
	  wave_wrf->temp[ix][iy] /= wsum;
	  wave_wrf->bg[ix][iy] /= wsum;
	  wave_wrf->pt[ix][iy] /= wsum;
	} else {
	  wave_wrf->temp[ix][iy] = GSL_NAN;
	  wave_wrf->bg[ix][iy] = GSL_NAN;
	  wave_wrf->pt[ix][iy] = GSL_NAN;
	}
      }

    /* Compute variance... */
    variance(wave_wrf, var_dh);

    /* Write WRF wave struct... */
    write_wave(argv[7], wave_wrf);
  }

  /* Free... */
  free(wave_airs);
  free(wave_wrf);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void background(
  double temp[NLAT][NLON],
  double pt[NLAT][NLON],
  int nlat,
  int nlon,
  int dlat,
  int dlon) {

  static double data[NLAT * NLAT];

  int ilon, ilat, ilon2, ilat2, n;

  /* Loop over grid points... */
  for (ilat = 0; ilat < nlat; ilat++)
    for (ilon = 0; ilon < nlon; ilon++) {

      /* Init... */
      n = 0;

      /* Average... */
      for (ilat2 = GSL_MAX(ilat - dlat, 0);
	   ilat2 <= GSL_MIN(ilat + dlat, nlat - 1); ilat2++)
	for (ilon2 = GSL_MAX(ilon - dlon, 0);
	     ilon2 <= GSL_MIN(ilon + dlon, nlon - 1); ilon2++)
	  if (gsl_finite(temp[ilat2][ilon2])) {
	    data[n] = temp[ilat2][ilon2];
	    n++;
	  }

      /* Set perturbation... */
      gsl_sort(data, 1, (size_t) n);
      pt[ilat][ilon] = temp[ilat][ilon]
	- gsl_stats_median_from_sorted_data(data, 1, (size_t) n);
    }
}
