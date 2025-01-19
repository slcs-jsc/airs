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
  Extract AIRS/Aqua orbit data.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;

  FILE *out;

  int i, track, xtrack;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG
      ("Give parameters: <orbit.tab> <airs_l1b_file> [ <airs_l1b_file2> ... ]");

  /* Create file... */
  printf("Write orbit data: %s\n", argv[1]);
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = satellite longitude [deg]\n"
	  "# $3 = satellite latitude [deg]\n"
	  "# $4 = footprint longitude [deg]\n"
	  "# $5 = footprint latitude [deg]\n");

  /* Loop over files... */
  for (i = 2; i < argc; i++) {

    /* Read AIRS data... */
    printf("Read AIRS Level-1B data file: %s\n", argv[i]);
    airs_rad_rdr(argv[i], &airs_rad_gran);

    /* Write data... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++) {
      fprintf(out, "\n");
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
	fprintf(out, "%.2f %g %g %g %g\n",
		airs_rad_gran.Time[track][xtrack] - 220838400,
		airs_rad_gran.sat_lon[track],
		airs_rad_gran.sat_lat[track],
		airs_rad_gran.Longitude[track][xtrack],
		airs_rad_gran.Latitude[track][xtrack]);
    }
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
