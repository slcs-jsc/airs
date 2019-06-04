#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum model dimensions. */
#define NLON 1441
#define NLAT 721
#define NZ 138

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static atm_t atm;
  static obs_t obs;

  static char pertname[LEN];

  static double lon[NLON], lat[NLAT], temp[NLAT][NLON], x0[3],
    x1[NLAT][NLON][3], w, wsum, rmax2 = 50. * 50., fwhm = 20., var_dh = 100.,
    btmean;

  static float z[NZ][NLAT][NLON], p[NZ][NLAT][NLON], t[NZ][NLAT][NLON],
    help[NZ * NLAT * NLON];

  static int id, itrack, ixtrack, ncid, dimid, varid,
    ilon, ilat, iz, nlon, nlat, nz;

  static size_t rs;

  pert_t *pert;

  wave_t *wave;

  FILE *out;

  /* ------------------------------------------------------------
     Get control parameters...
     ------------------------------------------------------------ */

  /* Check arguments... */
  if (argc < 8)
    ERRMSG("Give parameters: <ctl> <model.nc> <prof.tab> <map.tab> "
	   "<rad.tab> <pert.nc> <wave_airs.tab> <wave_model.tab>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);

  /* Set control parameters... */
  ctl.write_bbt = 1;

  /* Allocate... */
  ALLOC(pert, pert_t, 1);
  ALLOC(wave, wave_t, 1);

  /* ------------------------------------------------------------
     Read model data...
     ------------------------------------------------------------ */

  /* Open file... */
  printf("Read model data: %s\n", argv[2]);
  NC(nc_open(argv[2], NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "lev_2", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nz = (int) rs;
  if (nz > NZ)
    ERRMSG("Too many altitudes!");

  NC(nc_inq_dimid(ncid, "lat", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nlat = (int) rs;
  if (nlat > NLAT)
    ERRMSG("Too many latitudes!");

  NC(nc_inq_dimid(ncid, "lon", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nlon = (int) rs;
  if (nlon > NLON)
    ERRMSG("Too many longitudes!");

  /* Read latitudes... */
  NC(nc_inq_varid(ncid, "lat", &varid));
  NC(nc_get_var_double(ncid, varid, lat));

  /* Read longitudes... */
  NC(nc_inq_varid(ncid, "lon", &varid));
  NC(nc_get_var_double(ncid, varid, lon));

  /* Read temperature... */
  NC(nc_inq_varid(ncid, "t", &varid));
  NC(nc_get_var_float(ncid, varid, help));
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	t[iz][ilat][ilon] = help[(iz * nlat + ilat) * nlon + ilon];

  /* Read geopotential heights... */
  NC(nc_inq_varid(ncid, "gh", &varid));
  NC(nc_get_var_float(ncid, varid, help));
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	z[iz][ilat][ilon] =
	  (float) (help[(iz * nlat + ilat) * nlon + ilon] / 1e3);

  /* Calculate pressure... */
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	p[iz][ilat][ilon] = (float) (1013.25 * exp(-z[iz][ilat][ilon] / 7.0));

  /* Close file... */
  NC(nc_close(ncid));

  /* ------------------------------------------------------------
     Write model data to ASCII...
     ------------------------------------------------------------ */

  /* Write profile... */
  if (argv[3][0] != '-') {

    /* Create file... */
    printf("Write profile data: %s\n", argv[3]);
    if (!(out = fopen(argv[3], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = altitude index\n"
	    "# $2 = altitude [km]\n"
	    "# $3 = longitude [deg]\n"
	    "# $4 = latitude [deg]\n"
	    "# $5 = pressure [hPa]\n" "# $6 = temperature [K]\n\n");

    /* Write output... */
    for (iz = 0; iz < nz; iz++)
      fprintf(out, "%d %g %g %g %g %g\n", iz, z[iz][nlat / 2][nlon / 2],
	      lon[nlon / 2], lat[nlat / 2],
	      p[iz][nlat / 2][nlon / 2], t[iz][nlat / 2][nlon / 2]);

    /* Close file... */
    fclose(out);
  }

  /* Write map data... */
  if (argv[4][0] != '-') {

    /* Create file... */
    printf("Write map data: %s\n", argv[4]);
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
    for (ilon = 0; ilon < nlon; ilon++) {
      fprintf(out, "\n");
      for (ilat = 0; ilat < nlat; ilat++)
	fprintf(out, "%d %g %g %g %g %g\n", nz / 2, z[nz / 2][ilat][ilon],
		lon[ilon], lat[ilat],
		p[nz / 2][ilat][ilon], t[nz / 2][ilat][ilon]);
    }

    /* Close file... */
    fclose(out);
  }

  /* ------------------------------------------------------------
     Run forward model...
     ------------------------------------------------------------ */

  /* Set observation data... */
  obs.nr = 1;
  obs.obsz[0] = 700;

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
	    && z[iz][ilat][ilon] > 10 && z[iz][ilat][ilon] < 90) {
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
	    && z[iz][ilat][ilon] > 10 && z[iz][ilat][ilon] < 90) {
	  atm.p[atm.np] = p[iz][ilat][ilon];
	  atm.t[atm.np] = t[iz][ilat][ilon];
	  atm.np++;
	}

      /* Add top level... */
      atm.np++;

      /* Run forward model... */
      formod(&ctl, &atm, &obs);

      /* Get mean brightness temperature... */
      temp[ilat][ilon] = 0;
      for (id = 0; id < ctl.nd; id++)
	temp[ilat][ilon] += obs.rad[id][0] / ctl.nd;
    }
  }

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
	    "# $2 = latitude [deg]\n" "# $3 = brightness temperature [K]\n");

    /* Write output... */
    for (ilat = 0; ilat < nlat; ilat++) {
      fprintf(out, "\n");
      for (ilon = 0; ilon < nlon; ilon++)
	fprintf(out, "%g %g %g\n", lon[ilon], lat[ilat], temp[ilat][ilon]);
    }

    /* Close file... */
    fclose(out);
  }

  /* ------------------------------------------------------------
     Read AIRS radiance map and resample model data...
     ------------------------------------------------------------ */

  /* Check filename... */
  if (argv[6][0] != '-' && argv[7][0] != '-' && argv[8][0] != '-') {

    /* Read perturbation data... */
    read_pert(argv[6], pertname, pert);

    /* Convert to wave analysis struct... */
    pert2wave(pert, wave, 0, pert->ntrack - 1, 0, pert->nxtrack - 1);

    /* Estimate background... */
    background_poly(wave, 5, 0);

    /* Compute variance... */
    variance(wave, var_dh);

    /* Write observation wave struct... */
    write_wave(argv[7], wave);

    /* Get Cartesian coordinates for model grid... */
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	geo2cart(0, lon[ilon], lat[ilat], x1[ilat][ilon]);

    /* Loop over AIRS geolocations... */
    for (itrack = 0; itrack < pert->ntrack; itrack++)
      for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {

	/* Write info... */
	if (ixtrack == 0)
	  printf("  Average for track %d / %d ...\n",
		 itrack + 1, pert->ntrack);

	/* Init... */
	wsum = 0;
	btmean = 0;

	/* Average... */
	geo2cart(0, pert->lon[itrack][ixtrack],
		 pert->lat[itrack][ixtrack], x0);
	for (ilat = 0; ilat < nlat; ilat++)
	  for (ilon = 0; ilon < nlon; ilon++)
	    if (DIST2(x0, x1[ilat][ilon]) <= rmax2) {
	      w = exp(-DIST2(x0, x1[ilat][ilon]) /
		      (2. * gsl_pow_2(fwhm / 2.3548)));
	      btmean += w * temp[ilat][ilon];
	      wsum += w;
	    }

	/* Normalize... */
	if (wsum > 0)
	  pert->bt[itrack][ixtrack] = btmean / wsum;
	else
	  pert->bt[itrack][ixtrack] = GSL_NAN;
      }

    /* Convert to wave analysis struct... */
    pert2wave(pert, wave, 0, pert->ntrack - 1, 0, pert->nxtrack - 1);

    /* Estimate background... */
    background_poly(wave, 5, 0);

    /* Compute variance... */
    variance(wave, var_dh);

    /* Write model wave struct... */
    write_wave(argv[8], wave);
  }

  /* Free... */
  free(pert);
  free(wave);

  return EXIT_SUCCESS;
}
