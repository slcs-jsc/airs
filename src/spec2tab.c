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
  Extract radiance spectra.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;

  FILE *out;

  double dmin = 1e100, x0[3], x1[3];

  int ichan, track = -1, track2, xtrack = -1, xtrack2;

  /* Check arguments... */
  if (argc != 6)
    ERRMSG("Give parameters: <airs_l1b_file> "
	   "[index <track> <xtrack> | geo <lon> <lat>] <spec.tab>");

  /* Read AIRS data... */
  printf("Read AIRS Level-1B data file: %s\n", argv[1]);
  airs_rad_rdr(argv[1], &airs_rad_gran);

  /* Get indices... */
  if (argv[2][0] == 'i') {
    track = atoi(argv[3]);
    xtrack = atoi(argv[4]);
  }

  /* Find nearest footprint... */
  else {
    geo2cart(0, atof(argv[3]), atof(argv[4]), x0);
    for (track2 = 0; track2 < AIRS_RAD_GEOTRACK; track2++)
      for (xtrack2 = 0; xtrack2 < AIRS_RAD_GEOXTRACK; xtrack2++) {
	geo2cart(0, airs_rad_gran.Longitude[track2][xtrack2],
		 airs_rad_gran.Latitude[track2][xtrack2], x1);
	if (DIST2(x0, x1) < dmin) {
	  dmin = DIST2(x0, x1);
	  track = track2;
	  xtrack = xtrack2;
	}
      }
    if (dmin > 2500)
      ERRMSG("Geolocation not covered by granule!");
    printf("nearest footprint: lon= %g, lat= %g, track= %d, xtrack=%d\n",
	   airs_rad_gran.Longitude[track][xtrack],
	   airs_rad_gran.Latitude[track][xtrack], track, xtrack);
  }

  /* Check indices... */
  if (track < 0 || track >= AIRS_RAD_GEOTRACK)
    ERRMSG("Along-track index out of range!");
  if (xtrack < 0 || xtrack >= AIRS_RAD_GEOXTRACK)
    ERRMSG("Across-track index out of range!");

  /* Flag bad observations... */
  for (ichan = 0; ichan < AIRS_RAD_CHANNEL; ichan++)
    if ((airs_rad_gran.state[track][xtrack] != 0)
	|| (airs_rad_gran.ExcludedChans[ichan] > 2)
	|| (airs_rad_gran.CalChanSummary[ichan] & 8)
	|| (airs_rad_gran.CalChanSummary[ichan] & (32 + 64))
	|| (airs_rad_gran.CalFlag[track][ichan] & 16))
      airs_rad_gran.radiances[track][xtrack][ichan]
	= (float) sqrt(-1.0);

  /* Create file... */
  printf("Write spectrum: %s\n", argv[5]);
  if (!(out = fopen(argv[5], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = satellite longitude [deg]\n"
	  "# $3 = satellite latitude [deg]\n"
	  "# $4 = footprint longitude [deg]\n"
	  "# $5 = footprint latitude [deg]\n"
	  "# $6 = wavenumber [cm^-1]\n"
	  "# $7 = brightness temperature [K]\n"
	  "# $8 = radiance [W/(m^2 sr cm^-1)]\n\n");

  /* Write data... */
  for (ichan = 0; ichan < AIRS_RAD_CHANNEL; ichan++) {
    if (ichan > 0)
      if (fabs(airs_rad_gran.nominal_freq[ichan]
	       - airs_rad_gran.nominal_freq[ichan - 1]) > 1.2)
	fprintf(out, "\n");
    fprintf(out, "%.2f %g %g %g %g %g %g %g\n",
	    airs_rad_gran.Time[track][xtrack] - 220838400,
	    airs_rad_gran.sat_lon[track],
	    airs_rad_gran.sat_lat[track],
	    airs_rad_gran.Longitude[track][xtrack],
	    airs_rad_gran.Latitude[track][xtrack],
	    airs_rad_gran.nominal_freq[ichan],
	    BRIGHT(airs_rad_gran.radiances[track][xtrack][ichan] * 1e-3,
		   airs_rad_gran.nominal_freq[ichan]),
	    airs_rad_gran.radiances[track][xtrack][ichan] * 1e-3);
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
