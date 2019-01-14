#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of bands... */
#define NB 100

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  FILE *out;

  static airs_rad_gran_t airs_rad_gran;

  static double rad[NB];

  static int chan_min[NB], chan_max[NB], iarg, ib, ichan, n, nb, track,
    xtrack;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <out.tab> <l1b_file1> [<l1b_file2> ...]");

  /* Get control parameters... */
  nb = (int) scan_ctl(argc, argv, "NB", -1, "1", NULL);
  if (nb > NB)
    ERRMSG("Too many bands!");
  for (ib = 0; ib < nb; ib++) {
    chan_min[ib] = (int) scan_ctl(argc, argv, "CHAN_MIN", ib, "", NULL);
    if (chan_min[ib] < 0 || chan_min[ib] >= AIRS_RAD_CHANNEL)
      ERRMSG("Channel index out of range!");
    chan_max[ib] = (int) scan_ctl(argc, argv, "CHAN_MAX", ib, "", NULL);
    if (chan_max[ib] < 0 || chan_max[ib] >= AIRS_RAD_CHANNEL)
      ERRMSG("Channel index out of range!");
  }

  /* Create file... */
  printf("Write band data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Loop over HDF files... */
  for (iarg = 3; iarg < argc; iarg++) {

    /* Read AIRS data... */
    printf("Read AIRS Level-1B data file: %s\n", argv[iarg]);
    airs_rad_rdr(argv[iarg], &airs_rad_gran);

    /* Write header... */
    if (iarg == 3) {
      fprintf(out,
	      "# $1 = time [s]\n"
	      "# $2 = footprint longitude [deg]\n"
	      "# $3 = footprint latitude [deg]\n"
	      "# $4 = satellite altitude [km]\n"
	      "# $5 = satellite longitude [deg]\n"
	      "# $6 = satellite latitude [deg]\n");
      for (ib = 0; ib < nb; ib++)
	fprintf(out,
		"# $%d = BT(%.2f/cm...%.2f/cm) [K]\n",
		7 + ib, airs_rad_gran.nominal_freq[chan_min[ib]],
		airs_rad_gran.nominal_freq[chan_max[ib]]);
    }

    /* Flag bad observations... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
	for (ichan = 0; ichan < AIRS_RAD_CHANNEL; ichan++)
	  if ((airs_rad_gran.state[track][xtrack] != 0)
	      || (airs_rad_gran.ExcludedChans[ichan] > 2)
	      || (airs_rad_gran.CalChanSummary[ichan] & 8)
	      || (airs_rad_gran.CalChanSummary[ichan] & (32 + 64))
	      || (airs_rad_gran.CalFlag[track][ichan] & 16))
	    airs_rad_gran.radiances[track][xtrack][ichan] = GSL_NAN;

    /* Loop over scans... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++) {

      /* Write output... */
      fprintf(out, "\n");

      /* Loop over footprints... */
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {

	/* Write output... */
	fprintf(out, "%.2f %.4f %.4f %.3f %.4f %.4f",
		airs_rad_gran.Time[track][xtrack] - 220838400,
		airs_rad_gran.Longitude[track][xtrack],
		airs_rad_gran.Latitude[track][xtrack],
		airs_rad_gran.satheight[track],
		airs_rad_gran.sat_lon[track], airs_rad_gran.sat_lat[track]);

	/* Loop over bands... */
	for (ib = 0; ib < nb; ib++) {

	  /* Get mean radiance... */
	  n = 0;
	  rad[ib] = 0;
	  for (ichan = chan_min[ib]; ichan <= chan_max[ib]; ichan++)
	    if (gsl_finite(airs_rad_gran.radiances[track][xtrack][ichan])) {
	      rad[ib] += airs_rad_gran.radiances[track][xtrack][ichan];
	      n++;
	    }
	  if (n > 0)
	    rad[ib] /= n;
	  else
	    rad[ib] = GSL_NAN;

	  /* Convert to brightness temperature... */
	  rad[ib] = brightness(rad[ib] * 0.001,
			       0.5 *
			       (airs_rad_gran.nominal_freq[chan_min[ib]] +
				airs_rad_gran.nominal_freq[chan_max[ib]]));

	  /* Write output... */
	  fprintf(out, " %.3f", rad[ib]);
	}

	/* Write output... */
	fprintf(out, "\n");
      }
    }
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
