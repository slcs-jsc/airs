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
  Write retrieval data to ASCII file.
*/

#include "libairs.h"

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/* Replace dummy values by nan. */
#define CHECK(x) ((x)!=-9999 ? (x) : GSL_NAN)

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static airs_ret_gran_t airs_ret_gran;

  FILE *out;

  int lay, track, xtrack;

  /* Check arguments... */
  if (argc != 4)
    ERRMSG("Give parameters: <airs_l2_file> <layer> <airs.tab>");

  /* Get arguments... */
  lay = atoi(argv[2]);

  /* Read AIRS data... */
  printf("Read AIRS Level-2 data file: %s\n", argv[1]);
  airs_ret_rdr(argv[1], &airs_ret_gran);

  /* Create output file... */
  printf("Write ASCII file: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2  = altitude [km]\n"
	  "# $3  = longitude [deg]\n"
	  "# $4  = latitude [deg]\n"
	  "# $5  = pressure [hPa]\n"
	  "# $6  = temperature [K]\n"
	  "# $7  = H2O mass mixing ratio\n"
	  "# $8  = O3 volume mixing ratio\n"
	  "# $9  = CH4 volume mixing ratio\n"
	  "# $10 = CO volume mixing ratio\n");

  /* Write data to stdout... */
  for (track = 0; track < AIRS_RET_GEOTRACK; track++) {
    fprintf(out, "\n");
    for (xtrack = 0; xtrack < AIRS_RET_GEOXTRACK; xtrack++)
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g\n",
	      airs_ret_gran.Time[track][xtrack] - 220838400,
	      CHECK(airs_ret_gran.GP_Height[track][xtrack][lay]) / 1000,
	      CHECK(airs_ret_gran.Longitude[track][xtrack]),
	      CHECK(airs_ret_gran.Latitude[track][xtrack]),
	      CHECK(airs_ret_gran.pressStd[lay]),
	      CHECK(airs_ret_gran.TAirStd[track][xtrack][lay]),
	      CHECK(airs_ret_gran.H2OMMRStd[track][xtrack][lay]),
	      CHECK(airs_ret_gran.O3VMRStd[track][xtrack][lay]),
	      CHECK(airs_ret_gran.COVMRLevStd[track][xtrack][lay]),
	      CHECK(airs_ret_gran.CH4VMRLevStd[track][xtrack][lay]));
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
