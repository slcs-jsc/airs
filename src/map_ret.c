#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static ret_t ret;
  static wave_t wave;

  static double tbg[NDS], tabg[NDS], z0;

  FILE *out;

  char set[LEN];

  int asc, bg_poly_x, bg_poly_y, bg_smooth_x, bg_smooth_y, ids, ip, ix, iy;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <airs.nc> <map.tab>");

  /* Get control parameters... */
  scan_ctl(argc, argv, "SET", -1, "full", set);
  z0 = scan_ctl(argc, argv, "Z0", -1, "", NULL);
  bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "5", NULL);
  bg_poly_y = (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  bg_smooth_x = (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  bg_smooth_y = (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "0", NULL);

  /* Read AIRS data... */
  read_retr(argv[2], &ret);

  /* Get altitude index... */
  for (ip = 0; ip <= ret.np; ip++) {
    if (ip == ret.np)
      ERRMSG("Altitude level not found!");
    if (fabs(ret.z[0][ip] - z0) < 0.1)
      break;
  }

  /* Compute background... */
  ret2wave(&ret, &wave, 1, ip);
  background_poly(&wave, bg_poly_x, bg_poly_y);
  background_smooth(&wave, bg_smooth_x, bg_smooth_y);
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++)
      tbg[iy * 90 + ix] = wave.bg[ix][iy];
  ret2wave(&ret, &wave, 2, ip);
  background_poly(&wave, bg_poly_x, bg_poly_y);
  background_smooth(&wave, bg_smooth_x, bg_smooth_y);
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++)
      tabg[iy * 90 + ix] = wave.bg[ix][iy];

  /* Create output file... */
  printf("Write AIRS map data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2  = altitude [km]\n"
	  "# $3  = longitude [deg]\n"
	  "# $4  = latitude [deg]\n"
	  "# $5  = pressure [hPa]\n"
	  "# $6  = temperature (retrieved) [K]\n"
	  "# $7  = temperature (retrieved) perturbation [K]\n"
	  "# $8  = temperature (a priori) [K]\n"
	  "# $9  = temperature (a priori) perturbation [K]\n");
  fprintf(out,
	  "# $10 = temperature (total error) [K]\n"
	  "# $11 = temperature (noise error) [K]\n"
	  "# $12 = temperature (forward model error) [K]\n"
	  "# $13 = temperature (measurement content)\n"
	  "# $14 = temperature (resolution)\n" "# $15 = normalized chi^2\n");

  /* Write data... */
  for (ids = 0; ids < ret.nds; ids++) {

    /* Write new line... */
    if (ids % 90 == 0)
      fprintf(out, "\n");

    /* Check data... */
    if (ret.lon[ids][ip] < -180 || ret.lon[ids][ip] > 180
	|| ret.lat[ids][ip] < -90 || ret.lat[ids][ip] > 90
	|| ret.t[ids][ip] < 100 || ret.t[ids][ip] > 400)
      continue;

    /* Get ascending/descending flag... */
    asc = (ret.lat[ids > 90 ? ids : ids + 90][0]
	   > ret.lat[ids > 90 ? ids - 90 : ids][0]);

    /* Write data... */
    if (set[0] == 'f' || (set[0] == 'a' && asc) || (set[0] == 'd' && !asc))
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      ret.time[ids][ip], ret.z[ids][ip],
	      ret.lon[ids][ip], ret.lat[ids][ip],
	      ret.p[ids][ip], ret.t[ids][ip], ret.t[ids][ip] - tbg[ids],
	      ret.t_apr[ids][ip], ret.t_apr[ids][ip] - tabg[ids],
	      ret.t_tot[ids][ip], ret.t_noise[ids][ip], ret.t_fm[ids][ip],
	      ret.t_cont[ids][ip], ret.t_res[ids][ip], ret.chisq[ids]);
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
