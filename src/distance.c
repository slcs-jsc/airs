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
  Calculate distance between to geolocations.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  double x0[3], x1[3];

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <lon0> <lat0> <lon1> <lat1>");

  /* Read geolocations... */
  const double lon0 = atof(argv[1]);
  const double lat0 = atof(argv[2]);
  const double lon1 = atof(argv[3]);
  const double lat1 = atof(argv[4]);

  /* Write distance to stdout... */
  geo2cart(0, lon0, lat0, x0);
  geo2cart(0, lon1, lat1, x1);
  printf("%g\n", DIST(x0, x1));

  return EXIT_SUCCESS;
}
