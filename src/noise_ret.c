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
  Estimate noise based on retrieval data.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static ret_t ret;

  wave_t *wave, *wave2;

  FILE *out;

  double mu, mu2, nedt, nedt2;

  int ip;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <airs.nc> <noise.tab>");

  /* Allocate... */
  ALLOC(wave, wave_t, 1);
  ALLOC(wave2, wave_t, 1);

  /* Read AIRS data... */
  read_retr(argv[2], &ret);

  /* Create file... */
  printf("Write noise data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = altitude [km]\n"
	  "# $2 = longitude [deg]\n"
	  "# $3 = latitude [deg]\n"
	  "# $4 = mean temperature (retrieval) [K]\n"
	  "# $5 = noise estimate (retrieval) [K]\n"
	  "# $6 = mean temperature (a priori) [K]\n"
	  "# $7 = noise estimate (a priori) [K]\n\n");

  /* Loop over altitudes... */
  for (ip = 0; ip < ret.np; ip++) {

    /* Convert retrieval data to wave struct... */
    ret2wave(&ret, wave, 1, ip);
    ret2wave(&ret, wave2, 2, ip);

    /* Estimate noise... */
    noise(wave, &mu, &nedt);
    noise(wave2, &mu2, &nedt2);

    /* Estimate noise... */
    fprintf(out, "%g %g %g %g %g %g %g\n",
	    wave->z,
	    wave->lon[wave->nx / 2][wave->ny / 2],
	    wave->lat[wave->nx / 2][wave->ny / 2], mu, nedt, mu2, nedt2);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(wave);
  free(wave2);

  return EXIT_SUCCESS;
}
