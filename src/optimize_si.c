#include "libairs.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;

  static FILE *out;

  static double bt[AIRS_RAD_CHANNEL], bt2, dbt, lat0, lat1, lon0, lon1,
    mean[AIRS_RAD_CHANNEL][AIRS_RAD_CHANNEL],
    max[AIRS_RAD_CHANNEL][AIRS_RAD_CHANNEL],
    var[AIRS_RAD_CHANNEL][AIRS_RAD_CHANNEL];

  static int chan0, chan1, iarg, iavg, ichan, ichan2,
    n[AIRS_RAD_CHANNEL][AIRS_RAD_CHANNEL], navg, track, xtrack;

  /* Check arguments... */
  if (argc < 10)
    ERRMSG("Give parameters: <opt.tab> <chan0> <chan1>"
	   " <lon0> <lon1> <lat0> <lat1> <navg>"
	   " <l1b_file1> [<l1b_file2> ...]");

  /* Get parameters... */
  chan0 = GSL_MIN(GSL_MAX(atoi(argv[2]), 0), AIRS_RAD_CHANNEL - 1);
  chan1 = GSL_MIN(GSL_MAX(atoi(argv[3]), 0), AIRS_RAD_CHANNEL - 1);
  lon0 = atof(argv[4]);
  lon1 = atof(argv[5]);
  lat0 = atof(argv[6]);
  lat1 = atof(argv[7]);
  navg = atoi(argv[8]);

  /* Loop over HDF files... */
  for (iarg = 9; iarg < argc; iarg++) {

    /* Read AIRS data... */
    printf("Read AIRS Level-1B data file: %s\n", argv[iarg]);
    airs_rad_rdr(argv[iarg], &airs_rad_gran);

    /* Loop over footprints... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
	if (airs_rad_gran.Longitude[track][xtrack] >= lon0 &&
	    airs_rad_gran.Longitude[track][xtrack] <= lon1 &&
	    airs_rad_gran.Latitude[track][xtrack] >= lat0 &&
	    airs_rad_gran.Latitude[track][xtrack] <= lat1) {

	  /* Get brightness temperature... */
	  for (ichan = chan0; ichan <= chan1; ichan++)
	    if ((airs_rad_gran.state[track][xtrack] != 0)
		|| (airs_rad_gran.ExcludedChans[ichan] > 2)
		|| (airs_rad_gran.CalChanSummary[ichan] & 8)
		|| (airs_rad_gran.CalChanSummary[ichan] & (32 + 64))
		|| (airs_rad_gran.CalFlag[track][ichan] & 16))
	      bt[ichan] = GSL_NAN;
	    else
	      bt[ichan]
		= brightness(airs_rad_gran.radiances[track][xtrack][ichan]
			     * 0.001, airs_rad_gran.nominal_freq[ichan]);

	  /* Average channels... */
	  for (ichan = chan0; ichan <= chan1; ichan++) {
	    bt2 = 0;
	    for (iavg = 0; iavg < navg; iavg++)
	      bt2 += bt[ichan + iavg];
	    bt[ichan] = bt2 / navg;
	  }

	  /* Get statistics... */
	  for (ichan = chan0; ichan <= chan1; ichan++)
	    for (ichan2 = chan0; ichan2 <= chan1; ichan2++)
	      if (gsl_finite(bt[ichan]) && gsl_finite(bt[ichan2])) {

		/* Get brightness temperature difference... */
		dbt = (bt[ichan2] - bt[ichan]);
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
  for (ichan = chan0; ichan <= chan1; ichan++)
    for (ichan2 = chan0; ichan2 <= chan1; ichan2++) {
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
  for (ichan = chan0; ichan <= chan1; ichan++) {
    fprintf(out, "\n");
    for (ichan2 = chan0; ichan2 <= chan1; ichan2++)
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
