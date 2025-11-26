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
  Calculate solar zenith angle.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  /* Check arguments... */
  if (argc != 4)
    ERRMSG("Give parameters: <jsec> <lon> <lat>");

  /* Read arguments... */
  const double jsec = atof(argv[1]);
  const double lon = atof(argv[2]);
  const double lat = atof(argv[3]);

  /* Compute solar zenith angle... */
  printf("%g\n", RAD2DEG(acos(cos_sza(jsec, lon, lat))));
  
  return EXIT_SUCCESS;
}
