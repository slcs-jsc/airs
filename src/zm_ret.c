#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of latitudes. */
#define NLAT 180

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ret_t ret;
  static wave_t wave;

  static double apr_tm[NPG][NLAT], apr_var[NPG][NLAT], apr_noise[NPG][NLAT],
    ret_tm[NPG][NLAT], ret_var[NPG][NLAT], ret_noise[NPG][NLAT],
    ret_time[NPG][NLAT], mu, sig_apr, sig_ret, tbg[NDS], tabg[NDS];

  static int bg_poly_x, bg_poly_y, bg_smooth_x, bg_smooth_y,
    i, ids, ilat, ip, ix, iy, nlat, n[NPG][NLAT], ncid;

  FILE *out;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <zm.tab> <airs1.nc> [<airs2.nc> ...]");

  /* Get control parameters... */
  bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "5", NULL);
  bg_poly_y = (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  bg_smooth_x = (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  bg_smooth_y = (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "0", NULL);
  nlat = (int) scan_ctl(argc, argv, "NLAT", -1, "36", NULL);
  if (nlat > NLAT)
    ERRMSG("Too many latitudes!");

  /* Loop over files... */
  for (i = 3; i < argc; i++) {

    /* Read AIRS data... */
    if (nc_open(argv[i], NC_WRITE, &ncid) != NC_NOERR)
      continue;
    else
      nc_close(ncid);
    read_retr(argv[i], &ret);

    /* Loop over altitudes... */
    for (ip = 0; ip < ret.np; ip++) {

      /* Compute background... */
      ret2wave(&ret, &wave, 1, ip);
      background_poly(&wave, bg_poly_x, bg_poly_y);
      background_smooth(&wave, bg_smooth_x, bg_smooth_y);
      for (ix = 0; ix < wave.nx; ix++)
	for (iy = 0; iy < wave.ny; iy++)
	  tbg[iy * 90 + ix] = wave.bg[ix][iy];
      noise(&wave, &mu, &sig_ret);
      ret2wave(&ret, &wave, 2, ip);
      background_poly(&wave, bg_poly_x, bg_poly_y);
      background_smooth(&wave, bg_smooth_x, bg_smooth_y);
      for (ix = 0; ix < wave.nx; ix++)
	for (iy = 0; iy < wave.ny; iy++)
	  tabg[iy * 90 + ix] = wave.bg[ix][iy];
      noise(&wave, &mu, &sig_apr);

      /* Loop over data sets... */
      for (ids = 0; ids < ret.nds; ids++) {

	/* Check data... */
	if (ret.lon[ids][ip] < -180 || ret.lon[ids][ip] > 180
	    || ret.lat[ids][ip] < -90 || ret.lat[ids][ip] > 90
	    || ret.t[ids][ip] < 110 || ret.t[ids][ip] > 390
	    || !gsl_finite(ret.t[ids][ip]))
	  continue;

	/* Get latitude index... */
	ilat = (int) ((ret.lat[ids][ip] + 90.) / 180. * (double) nlat);
	if (ilat < 0 || ilat >= nlat)
	  continue;

	/* Get zonal mean... */
	if (gsl_finite(ret.t[ids][ip]) && gsl_finite(tbg[ids])) {
	  ret_time[ip][ilat] += ret.time[ids][ip];
	  ret_tm[ip][ilat] += ret.t[ids][ip];
	  ret_var[ip][ilat] += gsl_pow_2(ret.t[ids][ip] - tbg[ids]);
	  ret_noise[ip][ilat] += gsl_pow_2(sig_ret);
	  apr_tm[ip][ilat] += ret.t_apr[ids][ip];
	  apr_var[ip][ilat] += gsl_pow_2(ret.t_apr[ids][ip] - tabg[ids]);
	  apr_noise[ip][ilat] += gsl_pow_2(sig_apr);
	  n[ip][ilat]++;
	}
      }
    }
  }

  /* Create output file... */
  printf("Write AIRS zonal mean data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2  = altitude [km]\n"
	  "# $3  = latitude [deg]\n"
	  "# $4  = mean temperature (retrieved) [K]\n"
	  "# $5  = temperature variance (retrieved) [K^2]\n"
	  "# $6  = noise estimate (retrieved) [K^2]\n"
	  "# $7  = mean temperature (a priori) [K]\n"
	  "# $8  = temperature variance (a priori) [K^2]\n"
	  "# $9  = noise estimate (a priori) [K^2]\n"
	  "# $10 = number of data points\n");

  /* Loop over latitudes... */
  for (ilat = 0; ilat < nlat; ilat++) {

    /* Write empty line... */
    fprintf(out, "\n");

    /* Loop over altitudes... */
    for (ip = 0; ip < ret.np; ip++) {

      /* Write data... */
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %d\n",
	      ret_time[ip][ilat] / n[ip][ilat],
	      ret.z[0][ip], (ilat + 0.5) / nlat * 180. - 90.,
	      ret_tm[ip][ilat] / n[ip][ilat],
	      sqrt(ret_var[ip][ilat] / n[ip][ilat]),
	      sqrt(ret_noise[ip][ilat] / n[ip][ilat]),
	      apr_tm[ip][ilat] / n[ip][ilat],
	      sqrt(apr_var[ip][ilat] / n[ip][ilat]),
	      sqrt(apr_noise[ip][ilat] / n[ip][ilat]), n[ip][ilat]);
    }
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
