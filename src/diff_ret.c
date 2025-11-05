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
  Extract differences between retrieval data.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static retr_t ret, ret2;

  static FILE *out;

  static double mean[NPG], sigma[NPG], min[NPG], max[NPG],
    tt[NPG], lon[NPG], lat[NPG], temp[NPG], press[NPG];

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <airs.nc> <airs2.nc> <diff.tab>");

  /* Read AIRS data... */
  read_retr(argv[2], &ret);
  read_retr(argv[3], &ret2);

  /* Compute differences... */
  for (int ids = 0; ids < ret.nds; ids++)
    for (int ip = 0; ip < ret.np; ip++) {
      if (ret.time[ids][ip] != ret2.time[ids][ip] ||
	  ret.lon[ids][ip] != ret2.lon[ids][ip] ||
	  ret.lat[ids][ip] != ret2.lat[ids][ip])
	ERRMSG("Data files do not match!");
      tt[ip] += ret.time[ids][ip];
      lon[ip] += ret.lon[ids][ip];
      lat[ip] += ret.lat[ids][ip];
      press[ip] += ret.p[ids][ip];
      temp[ip] += ret.t[ids][ip];
      mean[ip] += ret2.t[ids][ip] - ret.t[ids][ip];
      sigma[ip] += gsl_pow_2(ret2.t[ids][ip] - ret.t[ids][ip]);
      min[ip] = GSL_MIN(min[ip], ret2.t[ids][ip] - ret.t[ids][ip]);
      max[ip] = GSL_MAX(max[ip], ret2.t[ids][ip] - ret.t[ids][ip]);
    }

  /* Create output file... */
  printf("Write retrieval differences data: %s\n", argv[4]);
  if (!(out = fopen(argv[4], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = pressure (set 1) [hPa]\n"
	  "# $6 = temperature (set 1) [K]\n"
	  "# $7 = temperature difference (mean, set 2 - set 1) [K]\n"
	  "# $8 = temperature difference (sigma, set 2 - set 1) [K]\n"
	  "# $9 = temperature difference (minimum, set 2 - set 1) [K]\n"
	  "# $10 = temperature difference (maximum, set 2 - set 1) [K]\n\n");

  /* Write output... */
  for (int ip = 0; ip < ret.np; ip++)
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g\n",
	    tt[ip] / ret.nds, ret.z[0][ip], lon[ip] / ret.nds,
	    lat[ip] / ret.nds, press[ip] / ret.nds, temp[ip] / ret.nds,
	    mean[ip] / ret.nds,
	    sqrt(sigma[ip] / ret.nds - gsl_pow_2(mean[ip] / ret.nds)),
	    min[ip], max[ip]);

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
