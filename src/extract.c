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
  Extract radiance data for retrievals.
*/

#include "libairs.h"

/* ------------------------------------------------------------
   Global variables...
   ------------------------------------------------------------ */

/* List of AIRS channels (don't change). */
int airs_chan[L1_NCHAN] = { 54, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
  2035, 2036, 2040, 2041, 2052, 2053, 2054, 2055,
  2067, 2075, 2076, 2077, 2078, 2079, 2080, 2081,
  2082, 2086, 2088, 2089, 2091, 2092, 2093
};

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Convert geopotential height to geometric altitude. */
double gph2z(
  double gph);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;
  static airs_ret_gran_t airs_ret_gran;

  static airs_l1_t l1;
  static airs_l2_t l2;

  /* Check arguments... */
  if (argc != 4)
    ERRMSG("Give parameters: <airs_l1_file> <airs_l2_file> <out.nc>");

  /* Check Level-1 filename... */
  if (argv[1][0] != '-') {

    /* Read data... */
    printf("Read AIRS Level-1 file: %s\n", argv[1]);
    airs_rad_rdr(argv[1], &airs_rad_gran);

    /* Flag bad data... */
    for (int track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (int xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
	for (int ichan = 0; ichan < L1_NCHAN; ichan++)
	  if ((airs_rad_gran.state[track][xtrack] != 0)
	      || (airs_rad_gran.ExcludedChans[airs_chan[ichan]] > 2)
	      || (airs_rad_gran.CalChanSummary[airs_chan[ichan]] & 8)
	      || (airs_rad_gran.CalChanSummary[airs_chan[ichan]] & (32 + 64))
	      || (airs_rad_gran.CalFlag[track][airs_chan[ichan]] & 16))
	    airs_rad_gran.radiances[track][xtrack][airs_chan[ichan]]
	      = GSL_NAN;

    /* Copy data to struct... */
    for (int track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (int xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
	l1.time[track][xtrack]
	  = airs_rad_gran.Time[track][xtrack] - 220838400.;
	l1.lon[track][xtrack]
	  = airs_rad_gran.Longitude[track][xtrack];
	l1.lat[track][xtrack]
	  = airs_rad_gran.Latitude[track][xtrack];
	l1.sat_z[track]
	  = airs_rad_gran.satheight[track];
	l1.sat_lon[track]
	  = airs_rad_gran.sat_lon[track];
	l1.sat_lat[track]
	  = airs_rad_gran.sat_lat[track];
	for (int ichan = 0; ichan < L1_NCHAN; ichan++) {
	  l1.nu[ichan]
	    = airs_rad_gran.nominal_freq[airs_chan[ichan]];
	  l1.rad[track][xtrack][ichan]
	    = airs_rad_gran.radiances[track][xtrack][airs_chan[ichan]] *
	    0.001f;
	}
      }

    /* Write netCDF file... */
    write_l1(argv[3], &l1);
  }

  /* Check Level-2 filename... */
  if (argv[2][0] != '-') {

    /* Read data... */
    printf("Read AIRS Level-2 file: %s\n", argv[2]);
    airs_ret_rdr(argv[2], &airs_ret_gran);

    /* Flag bad data... */
    for (int track = 0; track < AIRS_RET_GEOTRACK; track++)
      for (int xtrack = 0; xtrack < AIRS_RET_GEOXTRACK; xtrack++)
	for (int lay = 1; lay < AIRS_RET_STDPRESSURELAY; lay++)
	  if (airs_ret_gran.GP_Height[track][xtrack][lay] <= -9000.
	      || airs_ret_gran.TAirStd[track][xtrack][lay] <= -9000.) {
	    airs_ret_gran.GP_Height[track][xtrack][lay] = GSL_NAN;
	    airs_ret_gran.TAirStd[track][xtrack][lay] = GSL_NAN;
	  }

    /* Save data in struct... */
    for (int track = 0; track < AIRS_RET_GEOTRACK; track++)
      for (int xtrack = 0; xtrack < AIRS_RET_GEOXTRACK; xtrack++)
	for (int lay = 1; lay < AIRS_RET_STDPRESSURELAY; lay++) {
	  l2.time[track][xtrack]
	    = airs_ret_gran.Time[track][xtrack] - 220838400.;
	  l2.z[track][xtrack][lay - 1]
	    = airs_ret_gran.GP_Height[track][xtrack][lay] / 1000.;
	  l2.lon[track][xtrack]
	    = airs_ret_gran.Longitude[track][xtrack];
	  l2.lat[track][xtrack]
	    = airs_ret_gran.Latitude[track][xtrack];
	  l2.p[lay - 1]
	    = airs_ret_gran.pressStd[lay];
	  l2.t[track][xtrack][lay - 1]
	    = airs_ret_gran.TAirStd[track][xtrack][lay];
	}

    /* Convert geopotential heights to geometric heights... */
    for (int track = 0; track < L2_NTRACK; track++)
      for (int xtrack = 0; xtrack < L2_NXTRACK; xtrack++)
	for (int lay = 0; lay < L2_NLAY; lay++)
	  l2.z[track][xtrack][lay]
	    = gph2z(l2.z[track][xtrack][lay]);

    /* Write netCDF file... */
    write_l2(argv[3], &l2);
  }

  return EXIT_SUCCESS;
}

/*****************************************************************************/

double gph2z(
  double gph) {

  double a = 3.086e-3;

  return G0 / a - sqrt(gsl_pow_2(G0 / a) - 2 * G0 * gph / a);
}
