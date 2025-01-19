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
  Extract maps of retrieval data.
*/
#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static ret_t ret;
  static wave_t wave;

  static double tbg[NDS], tabg[NDS], z0;

  FILE *out;

  char set[LEN];

  int asc, bg_poly_x, bg_poly_y, bg_smooth_x, bg_smooth_y,
    ids, ip, ix, iy, npscan;

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
  npscan = (int) scan_ctl(argc, argv, "NPSCAN", -1, "90", NULL);

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
      tbg[iy * npscan + ix] = wave.bg[ix][iy];
  ret2wave(&ret, &wave, 2, ip);
  background_poly(&wave, bg_poly_x, bg_poly_y);
  background_smooth(&wave, bg_smooth_x, bg_smooth_y);
  for (ix = 0; ix < wave.nx; ix++)
    for (iy = 0; iy < wave.ny; iy++)
      tabg[iy * npscan + ix] = wave.bg[ix][iy];

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
    if (ids % npscan == 0)
      fprintf(out, "\n");

    /* Check data... */
    if (ret.lon[ids][ip] < -180 || ret.lon[ids][ip] > 180
	|| ret.lat[ids][ip] < -90 || ret.lat[ids][ip] > 90
	|| ret.t[ids][ip] < 100 || ret.t[ids][ip] > 400)
      continue;

    /* Get ascending/descending flag... */
    asc = (ret.lat[ids > npscan ? ids : ids + npscan][0]
	   > ret.lat[ids > npscan ? ids - npscan : ids][0]);

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
