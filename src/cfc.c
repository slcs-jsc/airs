#include "libairs.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Estimate noise. */
double get_noise(
  double bt,
  double dt250,
  double nu);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  FILE *out;

  static airs_rad_gran_t airs_rad_gran;

  static double ci, ci_err, ci_nedt = 0.35, cimax,
    f11_low, f11_low_err, f11_low_bt1, f11_low_bt1_nedt =
    0.35, f11_low_bt2, f11_low_bt2_nedt =
    0.32, f11_high, f11_high_err, f11_high_bt1, f11_high_bt1_nedt =
    0.34, f11_high_bt2, f11_high_bt2_nedt = 0.32;

  static int ichan, track, xtrack, iarg, f11_low_nu1 = 558, f11_low_nu2 =
    596, f11_high_nu1 = 624, f11_high_nu2 = 596, ci_nu = 558;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <out.tab> <l1b_file1> [<l1b_file2> ...]");

  /* Create file... */
  printf("Write CFC-11 emission data: %s\n", argv[1]);
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Loop over HDF files... */
  for (iarg = 2; iarg < argc; iarg++) {

    /* Read AIRS data... */
    printf("Read AIRS Level-1B data file: %s\n", argv[iarg]);
    airs_rad_rdr(argv[iarg], &airs_rad_gran);

    /* Write header... */
    if (iarg == 2) {
      fprintf(out,
	      "# $1  = time [s]\n"
	      "# $2  = footprint longitude [deg]\n"
	      "# $3  = footprint latitude [deg]\n"
	      "# $4  = satellite altitude [km]\n"
	      "# $5  = satellite longitude [deg]\n"
	      "# $6  = satellite latitude [deg]\n");
      fprintf(out,
	      "# $7  = cloud index, BT(%.2f/cm) [K]\n"
	      "# $8  = cloud index error [K]\n"
	      "# $9  = CFC-11 index (low wavenumbers),"
	      " BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $10 = CFC-11 index (low wavenumbers) error [K]\n"
	      "# $11 = CFC-11 index (high wavenumbers),"
	      " BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $12 = CFC-11 index (high wavenumbers) error [K]\n",
	      airs_rad_gran.nominal_freq[ci_nu],
	      airs_rad_gran.nominal_freq[f11_low_nu1],
	      airs_rad_gran.nominal_freq[f11_low_nu2],
	      airs_rad_gran.nominal_freq[f11_high_nu1],
	      airs_rad_gran.nominal_freq[f11_high_nu2]);
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

    /* Get maximum cloud index... */
    cimax = -999;
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
	ci = brightness(airs_rad_gran.radiances[track][xtrack][ci_nu] * 0.001,
			airs_rad_gran.nominal_freq[ci_nu]);
	if (ci > cimax)
	  cimax = ci;
      }

    /* Loop over scans... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++) {

      /* Write output... */
      fprintf(out, "\n");

      /* Loop over footprints... */
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {

	/* Skip daytime measurements... */
	if (sza(airs_rad_gran.Time[track][xtrack] - 220838400,
		airs_rad_gran.Longitude[track][xtrack],
		airs_rad_gran.Latitude[track][xtrack]) < 96.0)
	  continue;

	/* cloud index... */
	ci = brightness(airs_rad_gran.radiances[track][xtrack][ci_nu] * 0.001,
			airs_rad_gran.nominal_freq[ci_nu]);
	ci_err = get_noise(ci, ci_nedt, airs_rad_gran.nominal_freq[ci_nu]);

	/* Check cloud index... */
	if (ci < 0.95 * cimax || ci <= 270.)
	  continue;

	/* CFC-11 index (low wavenumbers)... */
	f11_low_bt1 =
	  brightness(airs_rad_gran.radiances[track][xtrack][f11_low_nu1] *
		     0.001, airs_rad_gran.nominal_freq[f11_low_nu1]);
	f11_low_bt2 =
	  brightness(airs_rad_gran.radiances[track][xtrack][f11_low_nu2] *
		     0.001, airs_rad_gran.nominal_freq[f11_low_nu2]);
	f11_low = f11_low_bt1 - f11_low_bt2;
	f11_low_err = sqrt(gsl_pow_2(get_noise(f11_low_bt1, f11_low_bt1_nedt,
					       airs_rad_gran.nominal_freq
					       [f11_low_nu1]))
			   +
			   gsl_pow_2(get_noise
				     (f11_low_bt2, f11_low_bt2_nedt,
				      airs_rad_gran.nominal_freq
				      [f11_low_nu2])));

	/* CFC-11 index (high wavenumbers)... */
	f11_high_bt1 =
	  brightness(airs_rad_gran.radiances[track][xtrack][f11_high_nu1] *
		     0.001, airs_rad_gran.nominal_freq[f11_high_nu1]);
	f11_high_bt2 =
	  brightness(airs_rad_gran.radiances[track][xtrack][f11_high_nu2] *
		     0.001, airs_rad_gran.nominal_freq[f11_high_nu2]);
	f11_high = f11_high_bt1 - f11_high_bt2;
	f11_high_err =
	  sqrt(gsl_pow_2
	       (get_noise
		(f11_high_bt1, f11_high_bt1_nedt,
		 airs_rad_gran.nominal_freq[f11_high_nu1]))
	       +
	       gsl_pow_2(get_noise
			 (f11_high_bt2, f11_high_bt2_nedt,
			  airs_rad_gran.nominal_freq[f11_high_nu2])));

	/* Write output... */
	fprintf(out,
		"%.2f %.4f %.4f %.3f %.4f %.4f %.2f %.2f %.2f %.2f %.2f %.2f\n",
		airs_rad_gran.Time[track][xtrack] - 220838400,
		airs_rad_gran.Longitude[track][xtrack],
		airs_rad_gran.Latitude[track][xtrack],
		airs_rad_gran.satheight[track],
		airs_rad_gran.sat_lon[track],
		airs_rad_gran.sat_lat[track],
		ci, ci_err, f11_low, f11_low_err, f11_high, f11_high_err);
      }
    }
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}

/************************************************************************/

double get_noise(
  double bt,
  double dt250,
  double nu) {

  double nesr;

  nesr = planck(250.0 + dt250, nu) - planck(250.0, nu);

  return brightness(planck(bt, nu) + nesr, nu) - bt;
}
