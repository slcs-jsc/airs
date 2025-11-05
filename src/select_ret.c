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
  Find retrieval results for given location.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static retr_t ret;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> [<airs1.nc> <airs2.c> ...]");

  /* Get control parameters... */
  const double lon0 =
    (int) scan_ctl(argc, argv, "SELECT_LON0", -1, "-180", NULL);
  const double lon1 =
    (int) scan_ctl(argc, argv, "SELECT_LON1", -1, "180", NULL);
  const double lat0 =
    (int) scan_ctl(argc, argv, "SELECT_LAT0", -1, "-90", NULL);
  const double lat1 =
    (int) scan_ctl(argc, argv, "SELECT_LAT1", -1, "90", NULL);

  /* Loop over retrieval files... */
  for (int i = 2; i < argc; i++) {

    /* Read AIRS data... */
    read_retr(argv[i], &ret);

    /* Check position... */
    if (ret.lon[ret.nds / 2][0] >= lon0
	&& ret.lon[ret.nds / 2][0] <= lon1
	&& ret.lat[ret.nds / 2][0] >= lat0 && ret.lat[ret.nds / 2][0] <= lat1)
      printf("select: %s %.2f %g %g\n", argv[i], ret.time[ret.nds / 2][0],
	     ret.lon[ret.nds / 2][0], ret.lat[ret.nds / 2][0]);
  }

  return EXIT_SUCCESS;
}
