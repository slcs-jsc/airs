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
  Optimize brightness temperature differences.
*/

#include "libairs.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;

  static FILE *out;

  static double bt[AIRS_RAD_CHANNEL],
    mean[AIRS_RAD_CHANNEL][AIRS_RAD_CHANNEL],
    max[AIRS_RAD_CHANNEL][AIRS_RAD_CHANNEL],
    var[AIRS_RAD_CHANNEL][AIRS_RAD_CHANNEL];

  static int n[AIRS_RAD_CHANNEL][AIRS_RAD_CHANNEL];

  /* Check arguments... */
  if (argc < 12)
    ERRMSG("Give parameters: <opt.tab> <sig_chan0> <sig_chan1>"
	   " <bg_chan0> <bg_chan1> <lon0> <lon1> <lat0> <lat1> <navg>"
	   " <l1b_file1> [<l1b_file2> ...]");

  /* Get parameters... */
  const int sig_chan0 =
    GSL_MIN(GSL_MAX(atoi(argv[2]), 0), AIRS_RAD_CHANNEL - 1);
  const int sig_chan1 =
    GSL_MIN(GSL_MAX(atoi(argv[3]), 0), AIRS_RAD_CHANNEL - 1);
  const int bg_chan0 =
    GSL_MIN(GSL_MAX(atoi(argv[4]), 0), AIRS_RAD_CHANNEL - 1);
  const int bg_chan1 =
    GSL_MIN(GSL_MAX(atoi(argv[5]), 0), AIRS_RAD_CHANNEL - 1);
  const double lon0 = atof(argv[6]);
  const double lon1 = atof(argv[7]);
  const double lat0 = atof(argv[8]);
  const double lat1 = atof(argv[9]);
  const int navg = atoi(argv[10]);

  /* Loop over HDF files... */
  for (int iarg = 11; iarg < argc; iarg++) {

    /* Read AIRS data... */
    printf("Read AIRS Level-1B data file: %s\n", argv[iarg]);
    airs_rad_rdr(argv[iarg], &airs_rad_gran);

    /* Loop over footprints... */
    for (int track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (int xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
	if (airs_rad_gran.Longitude[track][xtrack] >= lon0 &&
	    airs_rad_gran.Longitude[track][xtrack] <= lon1 &&
	    airs_rad_gran.Latitude[track][xtrack] >= lat0 &&
	    airs_rad_gran.Latitude[track][xtrack] <= lat1) {

	  /* Get brightness temperature... */
	  for (int ichan = 0; ichan < AIRS_RAD_CHANNEL; ichan++)
	    if ((airs_rad_gran.state[track][xtrack] != 0)
		|| (airs_rad_gran.ExcludedChans[ichan] > 2)
		|| (airs_rad_gran.CalChanSummary[ichan] & 8)
		|| (airs_rad_gran.CalChanSummary[ichan] & (32 + 64))
		|| (airs_rad_gran.CalFlag[track][ichan] & 16))
	      bt[ichan] = GSL_NAN;
	    else
	      bt[ichan]
		= BRIGHT(airs_rad_gran.radiances[track][xtrack][ichan]
			 * 0.001, airs_rad_gran.nominal_freq[ichan]);

	  /* Average channels... */
	  for (int ichan = 0; ichan < AIRS_RAD_CHANNEL - navg; ichan++) {
	    double bt2 = 0;
	    for (int iavg = 0; iavg < navg; iavg++)
	      bt2 += bt[ichan + iavg];
	    bt[ichan] = bt2 / navg;
	  }

	  /* Get statistics... */
	  for (int ichan = sig_chan0; ichan <= sig_chan1; ichan++)
	    for (int ichan2 = bg_chan0; ichan2 <= bg_chan1; ichan2++)
	      if (gsl_finite(bt[ichan]) && gsl_finite(bt[ichan2])) {

		/* Get brightness temperature difference... */
		const double dbt = (bt[ichan2] - bt[ichan]);
		if (fabs(dbt) > 100)
		  continue;

		/* Check filter... */
		if (n[ichan][ichan2] <= 0)
		  max[ichan][ichan2] = dbt;
		else
		  max[ichan][ichan2] = GSL_MAX(max[ichan][ichan2], dbt);
		mean[ichan][ichan2] += dbt;
		var[ichan][ichan2] += gsl_pow_2(dbt);
		n[ichan][ichan2]++;
	      }
	}
  }

  /* Normalize... */
  for (int ichan = sig_chan0; ichan <= sig_chan1; ichan++)
    for (int ichan2 = bg_chan0; ichan2 <= bg_chan1; ichan2++) {
      if (n[ichan][ichan2] > 0) {
	mean[ichan][ichan2] /= n[ichan][ichan2];
	var[ichan][ichan2] = sqrt(var[ichan][ichan2] / n[ichan][ichan2]
				  - gsl_pow_2(mean[ichan][ichan2]));
      } else
	mean[ichan][ichan2] = var[ichan][ichan2] = max[ichan][ichan2] =
	  GSL_NAN;
    }

  /* Write info... */
  printf("Write optimization data: %s\n", argv[1]);

  /* Create file... */
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = signal channel\n"
	  "# $2 = signal wavenumber [cm^-1]\n"
	  "# $3 = background channel\n"
	  "# $4 = background wavenumber [cm^-1]\n"
	  "# $5 = BTD(bg-sig) mean [K]\n"
	  "# $6 = BTD(bg-sig) standard deviation [K]\n"
	  "# $7 = BTD(bg-sig) maximum [K]\n"
	  "# $8 = effective SNR (= max/RMS)\n"
	  "# $9 = number of footprints\n");

  /* Write info... */
  for (int ichan = sig_chan0; ichan <= sig_chan1; ichan++) {
    fprintf(out, "\n");
    for (int ichan2 = bg_chan0; ichan2 <= bg_chan1; ichan2++)
      fprintf(out, "%d %.3f %d %.3f %g %g %g %g %d\n",
	      ichan, airs_rad_gran.nominal_freq[ichan],
	      ichan2, airs_rad_gran.nominal_freq[ichan2],
	      mean[ichan][ichan2], var[ichan][ichan2], max[ichan][ichan2],
	      max[ichan][ichan2] / sqrt(gsl_pow_2(var[ichan][ichan2])
					+ gsl_pow_2(mean[ichan][ichan2])),
	      n[ichan][ichan2]);
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
