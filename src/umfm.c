#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum UM dimensions. */
#define NLON 2310
#define NLAT 740
#define NZ 41

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Estimate background... */
void background(
  double temp[NLAT][NLON],
  double pt[NLAT][NLON],
  int nlat,
  int nlon,
  int bg_poly_x,
  int bg_smooth_y);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static atm_t atm;
  static obs_t obs;

  static double z[NZ], p[NZ][NLAT][NLON], t[NZ][NLAT][NLON],
    lon[NLON], lat[NLAT], temp[NLAT][NLON], pt[NLAT][NLON],
    x0[3], x1[NLAT][NLON][3], wsum, rmax2 = 10. * 10., var_dh;

  static int bg_poly_x, bg_smooth_y, id, ix, iy, oit, oiz,
    ncid, dimid, varid, ilon, ilon2, ilat, ilat2, iz, ncrop, nlon, nlat, nz;

  static size_t start[10], count[10], rs;

  wave_t *wave_airs, *wave_um;

  FILE *out;

  /* ------------------------------------------------------------
     Get control parameters...
     ------------------------------------------------------------ */

  /* Check arguments... */
  if (argc < 10)
    ERRMSG("Give parameters: <ctl> <ump.nc> <umtheta.nc> <it> "
	   "<wave_airs.tab> <out_um.tab> <iz> <out_rad.tab> <wave_um.tab>");

  /* Get arguments... */
  oit = atoi(argv[4]);
  oiz = atoi(argv[7]);

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Set control parameters... */
  ctl.write_bbt = 1;

  /* Get control parameters... */
  bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "5", NULL);
  bg_smooth_y = (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "7", NULL);
  ncrop = (int) scan_ctl(argc, argv, "NCROP", -1, "10", NULL);
  var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "0", NULL);

  /* Allocate... */
  ALLOC(wave_airs, wave_t, 1);
  ALLOC(wave_um, wave_t, 1);

  /* ------------------------------------------------------------
     Read UM data...
     ------------------------------------------------------------ */

  /* Read pressure file... */
  printf("Read UM pressure data: %s\n", argv[2]);
  NC(nc_open(argv[2], NC_NOWRITE, &ncid));

  /* Read longitudes... */
  NC(nc_inq_dimid(ncid, "longitude", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nlon = (int) rs;
  if (nlon >= NLON)
    ERRMSG("Too many longitudes!");
  NC(nc_inq_varid(ncid, "longitude", &varid));
  NC(nc_get_var_double(ncid, varid, lon));

  /* Read latitudes... */
  NC(nc_inq_dimid(ncid, "latitude", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nlat = (int) rs;
  if (nlat >= NLAT)
    ERRMSG("Too many latitudes!");
  NC(nc_inq_varid(ncid, "latitude", &varid));
  NC(nc_get_var_double(ncid, varid, lat));

  /* Read heights... */
  NC(nc_inq_dimid(ncid, "ht", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  nz = (int) rs;
  if (nz >= NZ)
    ERRMSG("Too many heights!");
  NC(nc_inq_varid(ncid, "ht", &varid));
  NC(nc_get_var_double(ncid, varid, z));

  /* Read pressure... */
  NC(nc_inq_varid(ncid, "p", &varid));
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

  /* Close file... */
  NC(nc_close(ncid));

  /* Read theta file... */
  printf("Read UM theta data: %s\n", argv[3]);
  NC(nc_open(argv[3], NC_NOWRITE, &ncid));

  /* Read theta... */
  NC(nc_inq_varid(ncid, "theta", &varid));
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

  /* Close file... */
  NC(nc_close(ncid));

  /* ------------------------------------------------------------
     Convert UM data...
     ------------------------------------------------------------ */

  /* Modify longitudes... */
  for (ilon = 0; ilon < nlon; ilon++)
    if (lon[ilon] > 180)
      lon[ilon] -= 360;

  /* Scale heights... */
  for (iz = 0; iz < nz; iz++)
    z[iz] /= 1e3;

  /* Scale pressure and theta... */
  for (iz = 0; iz < nz; iz++)
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	if (p[iz][ilat][ilon] <= 0 || p[iz][ilat][ilon] >= 1000000 ||
	    t[iz][ilat][ilon] <= 0 || t[iz][ilat][ilon] >= 10000) {
	  p[iz][ilat][ilon] = GSL_NAN;
	  t[iz][ilat][ilon] = GSL_NAN;
	} else {
	  p[iz][ilat][ilon] /= 1e2;
	  t[iz][ilat][ilon] /= pow(1e3 / p[iz][ilat][ilon], 0.286);
	}

  /* ------------------------------------------------------------
     Write UM data to ASCII...
     ------------------------------------------------------------ */

  /* Check filename... */
  if (argv[6][0] != '-') {

    /* Check height level... */
    if (oiz < 0 || oiz >= nz)
      ERRMSG("Height index out of range!");

    /* Create file... */
    printf("Write UM data: %s\n", argv[6]);
    if (!(out = fopen(argv[6], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = altitude [km]\n"
	    "# $2 = longitude [deg]\n"
	    "# $3 = latitude [deg]\n"
	    "# $4 = pressure [hPa]\n" "# $5 = temperature [K]\n");

    /* Write output... */
    for (ilon = 0; ilon < nlon; ilon++) {
      fprintf(out, "\n");
      for (ilat = 0; ilat < nlat; ilat++)
	fprintf(out, "%g %g %g %g %g\n", z[oiz], lon[ilon], lat[ilat],
		p[oiz][ilat][ilon], t[oiz][ilat][ilon]);
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

      /* Set atmospheric data... */
      atm.np = 0;
      for (iz = 0; iz < nz; iz++)
	if (gsl_finite(p[iz][ilat][ilon]) && gsl_finite(t[iz][ilat][ilon])
	    && p[iz][ilat][ilon] > 0 && p[iz][ilat][ilon] < 1200
	    && t[iz][ilat][ilon] > 100 && t[iz][ilat][ilon] < 400) {
	  atm.z[atm.np] = z[iz];
	  if ((++atm.np) >= NP)
	    ERRMSG("Too many altitudes!");
	}
      climatology(&ctl, &atm);
      atm.np = 0;
      for (iz = 0; iz < nz; iz++)
	if (gsl_finite(p[iz][ilat][ilon]) && gsl_finite(t[iz][ilat][ilon])
	    && p[iz][ilat][ilon] > 0 && p[iz][ilat][ilon] < 1200
	    && t[iz][ilat][ilon] > 100 && t[iz][ilat][ilon] < 400) {
	  atm.p[atm.np] = p[iz][ilat][ilon];
	  atm.t[atm.np] = t[iz][ilat][ilon];
	  atm.np++;
	}

      /* Check number of altitudes... */
      if (atm.np < 20) {
	temp[ilat][ilon] = GSL_NAN;
	continue;
      }

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
  for (ilat = 0; ilat < nlat; ilat++) {
    for (ilon = 0; ilon < nlon; ilon++)
      if (gsl_finite(temp[ilat][ilon])) {
	for (ilon2 = ilon; ilon2 <= GSL_MIN(ilon + ncrop, nlon - 1); ilon2++)
	  temp[ilat][ilon2] = GSL_NAN;
	break;
      }
    for (ilon = nlon - 1; ilon >= 0; ilon--)
      if (gsl_finite(temp[ilat][ilon])) {
	for (ilon2 = ilon; ilon2 >= GSL_MAX(ilon - ncrop, 0); ilon2--)
	  temp[ilat][ilon2] = GSL_NAN;
	break;
      }
  }
  for (ilon = 0; ilon < nlon; ilon++) {
    for (ilat = 0; ilat < nlat; ilat++)
      if (gsl_finite(temp[ilat][ilon])) {
	for (ilat2 = ilat; ilat2 <= GSL_MIN(ilat + ncrop, nlat - 1); ilat2++)
	  temp[ilat2][ilon] = GSL_NAN;
	break;
      }
    for (ilat = nlat - 1; ilat >= 0; ilat--)
      if (gsl_finite(temp[ilat][ilon])) {
	for (ilat2 = ilat; ilat2 >= GSL_MAX(ilat - ncrop, 0); ilat2--)
	  temp[ilat2][ilon] = GSL_NAN;
	break;
      }
  }

  /* Get perturbations... */
  background(temp, pt, nlat, nlon, bg_poly_x, bg_smooth_y);

  /* ------------------------------------------------------------
     Save forward model output...
     ------------------------------------------------------------ */

  /* Check filename... */
  if (argv[8][0] != '-') {

    /* Create file... */
    printf("Write radiance data: %s\n", argv[8]);
    if (!(out = fopen(argv[8], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = longitude [deg]\n"
	    "# $2 = latitude [deg]\n"
	    "# $3 = UM brightness temperature [K]\n"
	    "# $4 = UM brightness temperature perturbation [K]\n");

    /* Write output... */
    for (ilat = 0; ilat < nlat; ilat++) {
      fprintf(out, "\n");
      for (ilon = 0; ilon < nlon; ilon++)
	fprintf(out, "%g %g %g %g\n", lon[ilon], lat[ilat],
		temp[ilat][ilon], pt[ilat][ilon]);
    }

    /* Close file... */
    fclose(out);
  }

  /* ------------------------------------------------------------
     Read AIRS radiance map and resample model data...
     ------------------------------------------------------------ */

  /* Check filename... */
  if (argv[5][0] != '-') {

    /* Read AIRS wave file... */
    read_wave(argv[5], wave_airs);
    memcpy(wave_um, wave_airs, sizeof(wave_t));

    /* Get Cartesian coordinates for model grid... */
    for (ilat = 0; ilat < nlat; ilat++)
      for (ilon = 0; ilon < nlon; ilon++)
	geo2cart(0, lon[ilon], lat[ilat], x1[ilat][ilon]);

    /* Loop over AIRS geolocations... */
    for (ix = 0; ix < wave_airs->nx; ix++)
      for (iy = 0; iy < wave_airs->ny; iy++) {

	/* Write info... */
	if (iy == 0)
	  printf("  Average for xtrack %d / %d ...\n", ix + 1, wave_airs->nx);

	/* Init... */
	wsum = 0;
	wave_um->temp[ix][iy] = 0;
	wave_um->bg[ix][iy] = 0;
	wave_um->pt[ix][iy] = 0;
	wave_um->var[ix][iy] = 0;

	/* Average... */
	geo2cart(0, wave_airs->lon[ix][iy], wave_airs->lat[ix][iy], x0);
	for (ilat = 0; ilat < nlat; ilat++)
	  for (ilon = 0; ilon < nlon; ilon++)
	    if (DIST2(x0, x1[ilat][ilon]) <= rmax2) {
	      wave_um->temp[ix][iy] += temp[ilat][ilon];
	      wave_um->bg[ix][iy] += temp[ilat][ilon] - pt[ilat][ilon];
	      wave_um->pt[ix][iy] += pt[ilat][ilon];
	      wsum++;
	    }

	/* Normalize... */
	wave_um->temp[ix][iy] /= wsum;
	wave_um->bg[ix][iy] /= wsum;
	wave_um->pt[ix][iy] /= wsum;
      }

    /* Compute variance... */
    variance(wave_um, var_dh);

    /* Write UM wave struct... */
    write_wave(argv[9], wave_um);
  }

  /* Free... */
  free(wave_airs);
  free(wave_um);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void background(
  double temp[NLAT][NLON],
  double pt[NLAT][NLON],
  int nlat,
  int nlon,
  int bg_poly_x,
  int bg_smooth_y) {

  static double bg[NLAT][NLON];

  gsl_multifit_linear_workspace *work;
  gsl_matrix *cov, *X;
  gsl_vector *c, *x, *y;

  double chisq, bsum, wsum;

  int ilon, ilat, dlat;

  size_t dim, i, i2, n;

  /* Compute background... */
  for (ilat = 0; ilat < nlat; ilat++) {

    /* Get number of points... */
    n = 0;
    for (ilon = 0; ilon < nlon; ilon++) {
      bg[ilat][ilon] = GSL_NAN;
      if (gsl_finite(temp[ilat][ilon]))
	n++;
    }
    if (n < 10)
      continue;

    /* Allocate... */
    dim = (size_t) bg_poly_x;
    work = gsl_multifit_linear_alloc(n, dim);
    cov = gsl_matrix_alloc(dim, dim);
    X = gsl_matrix_alloc(n, dim);
    c = gsl_vector_alloc(dim);
    x = gsl_vector_alloc(n);
    y = gsl_vector_alloc(n);

    /* Fit polynomial... */
    i = 0;
    for (ilon = 0; ilon < nlon; ilon++)
      if (gsl_finite(temp[ilat][ilon])) {
	gsl_vector_set(x, i, (double) i);
	gsl_vector_set(y, i, temp[ilat][ilon]);
	for (i2 = 0; i2 < dim; i2++)
	  gsl_matrix_set(X, i, i2, pow(gsl_vector_get(x, i), (double) i2));
	i++;
      }
    gsl_multifit_linear(X, y, c, cov, &chisq, work);
    i = 0;
    for (ilon = 0; ilon < nlon; ilon++)
      if (gsl_finite(temp[ilat][ilon])) {
	bg[ilat][ilon] =
	  gsl_poly_eval(c->data, (int) dim, gsl_vector_get(x, i));
	i++;
      }

    /* Free... */
    gsl_multifit_linear_free(work);
    gsl_matrix_free(cov);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_vector_free(x);
    gsl_vector_free(y);
  }

  /* Smooth background and calculate perturbations... */
  for (ilon = 0; ilon < nlon; ilon++)
    for (ilat = 0; ilat < nlat; ilat++) {

      /* Smooth background... */
      bsum = wsum = 0;
      for (dlat = -bg_smooth_y; dlat <= bg_smooth_y; dlat++)
	if (ilat + dlat >= 0 && ilat + dlat < nlat) {
	  bsum += bg[ilat + dlat][ilon];
	  wsum++;
	}

      /* Compute perturbations... */
      pt[ilat][ilon] = temp[ilat][ilon] - bsum / wsum;
    }
}
