#include "libairs.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Estimate noise. */
double get_noise(
  double bt,
  double dt250,
  double nu);

/* Estimate SO2 column density. */
void get_so2_column(
  double si,
  double dsi,
  double t,
  double lat,
  int set,
  double *scd,
  double *err);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  FILE *out;

  static airs_rad_gran_t airs_rad_gran;

  static double ci, ci_err, ci_nedt = 0.0783,
    ai, ai_err, ai_bt1, ai_bt1_nedt = 0.3155, ai_bt2, ai_bt2_nedt = 0.1177,
    si_low, si_low_err, si_low_bt1, si_low_bt1_nedt = 0.1064,
    si_low_bt2, si_low_bt2_nedt = 0.0909,
    si_high, si_high_err, si_high_bt1, si_high_bt1_nedt = 0.1064,
    si_high_bt2, si_high_bt2_nedt = 0.0786,
    scd_low, scd_low_err, scd_high, scd_high_err, scd, scd_err;

  static int ichan, track, xtrack, iarg, ai_nu1 = 559, ai_nu2 = 901, ci_nu =
    1290, si_low_nu1 = 1591, si_low_nu2 = 1526, si_high_nu1 =
    1591, si_high_nu2 = 1550;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <out.tab> <l1b_file1> [<l1b_file2> ...]");

  /* Create file... */
  printf("Write volcanic emission data: %s\n", argv[1]);
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
	      "# $4  = cloud index, BT(%.2f/cm) [K]\n"
	      "# $5  = cloud index error [K]\n"
	      "# $6  = ash index, BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $7  = ash index error [K]\n",
	      airs_rad_gran.nominal_freq[ci_nu],
	      airs_rad_gran.nominal_freq[ai_nu1],
	      airs_rad_gran.nominal_freq[ai_nu2]);
      fprintf(out,
	      "# $8  = SO2 index (low), BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $9  = SO2 index (low) error [K]\n"
	      "# $10 = SO2 index (high), BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $11 = SO2 index (high) error [K]\n"
	      "# $12 = SO2 column density estimate [DU]\n"
	      "# $13 = SO2 column density error [DU]\n",
	      airs_rad_gran.nominal_freq[si_low_nu1],
	      airs_rad_gran.nominal_freq[si_low_nu2],
	      airs_rad_gran.nominal_freq[si_high_nu1],
	      airs_rad_gran.nominal_freq[si_high_nu2]);
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

	/* cloud index... */
	ci = brightness(airs_rad_gran.radiances[track][xtrack][ci_nu] * 0.001,
			airs_rad_gran.nominal_freq[ci_nu]);
	ci_err = get_noise(ci, ci_nedt, airs_rad_gran.nominal_freq[ci_nu]);

	/* ash index... */
	ai_bt1 =
	  brightness(airs_rad_gran.radiances[track][xtrack][ai_nu1] * 0.001,
		     airs_rad_gran.nominal_freq[ai_nu1]);
	ai_bt2 =
	  brightness(airs_rad_gran.radiances[track][xtrack][ai_nu2] * 0.001,
		     airs_rad_gran.nominal_freq[ai_nu2]);
	ai = ai_bt1 - ai_bt2;
	ai_err = sqrt(gsl_pow_2(get_noise(ai_bt1, ai_bt1_nedt,
					  airs_rad_gran.nominal_freq[ai_nu1]))
		      + gsl_pow_2(get_noise(ai_bt2, ai_bt2_nedt,
					    airs_rad_gran.nominal_freq
					    [ai_nu2])));

	/* SO2 index (low concentrations)... */
	si_low_bt1 =
	  brightness(airs_rad_gran.radiances[track][xtrack][si_low_nu1] *
		     0.001, airs_rad_gran.nominal_freq[si_low_nu1]);
	si_low_bt2 =
	  brightness(airs_rad_gran.radiances[track][xtrack][si_low_nu2] *
		     0.001, airs_rad_gran.nominal_freq[si_low_nu2]);
	si_low = si_low_bt1 - si_low_bt2;
	si_low_err = sqrt(gsl_pow_2(get_noise(si_low_bt1, si_low_bt1_nedt,
					      airs_rad_gran.nominal_freq
					      [si_low_nu1]))
			  +
			  gsl_pow_2(get_noise
				    (si_low_bt2, si_low_bt2_nedt,
				     airs_rad_gran.nominal_freq
				     [si_low_nu2])));

	/* SO2 index (high concentrations)... */
	si_high_bt1 =
	  brightness(airs_rad_gran.radiances[track][xtrack][si_high_nu1] *
		     0.001, airs_rad_gran.nominal_freq[si_high_nu1]);
	si_high_bt2 =
	  brightness(airs_rad_gran.radiances[track][xtrack][si_high_nu2] *
		     0.001, airs_rad_gran.nominal_freq[si_high_nu2]);
	si_high = si_high_bt1 - si_high_bt2;
	si_high_err = sqrt(gsl_pow_2(get_noise(si_high_bt1, si_high_bt1_nedt,
					       airs_rad_gran.nominal_freq
					       [si_high_nu1]))
			   +
			   gsl_pow_2(get_noise
				     (si_high_bt2, si_high_bt2_nedt,
				      airs_rad_gran.nominal_freq
				      [si_high_nu2])));

	/* SO2 column density (low concentrations)... */
	get_so2_column(si_low, si_low_err,
		       airs_rad_gran.Time[track][xtrack] - 220838400,
		       airs_rad_gran.Latitude[track][xtrack],
		       1, &scd_low, &scd_low_err);

	/* SO2 column density (high concentrations)... */
	get_so2_column(si_high, si_high_err,
		       airs_rad_gran.Time[track][xtrack] - 220838400,
		       airs_rad_gran.Latitude[track][xtrack],
		       2, &scd_high, &scd_high_err);

	/* Get optimal estimate... */
	scd =
	  (scd_low * gsl_pow_2(scd_high_err) +
	   scd_high * gsl_pow_2(scd_low_err))
	  / (gsl_pow_2(scd_low_err) + gsl_pow_2(scd_high_err));
	scd_err =
	  1 / sqrt(1 / gsl_pow_2(scd_low_err) + 1 / gsl_pow_2(scd_high_err));

	/* Write output... */
	fprintf(out,
		"%.2f %.4f %.4f %.2f %.2f %.2f %.2f "
		"%.2f %.2f %.2f %.2f %.1f %.1f\n",
		airs_rad_gran.Time[track][xtrack] - 220838400,
		airs_rad_gran.Longitude[track][xtrack],
		airs_rad_gran.Latitude[track][xtrack],
		ci, ci_err, GSL_MAX(ai, 0.0), ai_err,
		GSL_MAX(si_low, 0.0), si_low_err,
		GSL_MAX(si_high, 0.0), si_high_err, scd, scd_err);
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

/************************************************************************/

void get_so2_column(
  double si,
  double dsi,
  double t,
  double lat,
  int set,
  double *scd,
  double *err) {

  static double
    si_low[53] = { -0.377, -0.361, -0.342, -0.318, -0.291, -0.257, -0.217,
    -0.169, -0.112, -0.043, 0.039, 0.138, 0.256, 0.397,
    0.565, 0.766, 1.005, 1.29, 1.629, 2.03, 2.505, 3.065,
    3.725, 4.496, 5.398, 6.44, 7.644, 9.019, 10.574, 12.329,
    14.254, 16.378, 18.638, 21.039, 23.504, 25.989, 28.413,
    30.71, 32.786, 34.622, 36.118, 37.338, 38.216, 38.865,
    39.43, 39.886, 39.741, 39.86, 39.821, 39.832, 39.776,
    39.649, 39.659
  };

  static double
    scd_low[53] = { 0.205917, 0.232053, 0.263417, 0.301053, 0.346217,
    0.400413, 0.465446, 0.543491, 0.637141, 0.749524,
    0.884383, 1.04621, 1.24041, 1.47344, 1.75308,
    2.08865, 2.49133, 2.97455, 3.55441, 4.25026, 5.08524,
    6.08725, 7.28967, 8.73257, 10.464, 12.5418, 15.035,
    18.0271, 21.6174, 25.9259, 31.0959, 37.3, 44.745,
    53.6792, 64.4, 77.2647, 92.7026, 111.228, 133.458,
    160.135, 192.147, 230.562, 276.659, 331.977, 398.357,
    478.011, 1189.33, 1427.18, 2959.33, 3551.19, 5113.68,
    8836.36, 10603.6
  };

  static double
    si_high[60] = { -4.203, -4.199, -4.195, -4.19, -4.184, -4.177, -4.168,
    -4.158, -4.145, -4.13, -4.112, -4.091, -4.065, -4.034,
    -3.996, -3.952, -3.898, -3.834, -3.758, -3.666, -3.557,
    -3.426, -3.27, -3.084, -2.863, -2.599, -2.287, -1.918,
    -1.481, -0.966, -0.363, 0.343, 1.16, 2.107, 3.19, 4.421,
    5.811, 7.35, 9.049, 10.887, 12.852, 14.93, 17.065,
    19.269, 21.482, 23.711, 25.909, 28.064, 30.136, 32.094,
    33.877, 35.466, 36.773, 37.835, 38.59, 39.314, 39.866,
    39.826, 39.737, 39.791
  };

  static double
    scd_high[60] = { 0.205917, 0.232053, 0.263417, 0.301053, 0.346217,
    0.400413, 0.465446, 0.543491, 0.637141, 0.749524,
    0.884383, 1.04621, 1.24041, 1.47344, 1.75308, 2.08865,
    2.49133, 2.97455, 3.55441, 4.25026, 5.08524, 6.08725,
    7.28967, 8.73257, 10.464, 12.5418, 15.035, 18.0271,
    21.6174, 25.9259, 31.0959, 37.3, 44.745, 53.6792,
    64.4, 77.2647, 92.7026, 111.228, 133.458, 160.135,
    192.147, 230.562, 276.659, 331.977, 398.357, 478.011,
    573.599, 688.305, 825.952, 991.126, 1189.33, 1427.18,
    1712.61, 2055.12, 2466.13, 2959.33, 3551.19, 5113.68,
    7363.64, 10603.6
  };

  double *sia, *scda, scdm, scdp, s1, w_eqn, w_midl, w_psum, w_pwin;

  int i, *n, n_low = 53, n_high = 60;

  /* Set data set... */
  if (set == 1) {
    sia = &si_low[0];
    scda = &scd_low[0];
    n = &n_low;
  } else if (set == 2) {
    sia = &si_high[0];
    scda = &scd_high[0];
    n = &n_high;
  } else
    ERRMSG("Coding error!");

  /* Get weighting factors... */
  if (fabs(lat) <= 45) {
    w_eqn = LIN(0.0, 1.0, 45.0, 0.0, fabs(lat));
    w_midl = 1 - w_eqn;
    w_psum = 0;
    w_pwin = 0;
  } else {
    w_eqn = 0;
    w_midl = LIN(45.0, 1.0, 90.0, 0.0, fabs(lat));
    if (lat > 0) {
      w_psum = 0.5 * (1 - cos(2 * M_PI * t / (86400.0 * 365.25)));
      w_pwin = 1 - w_psum;
    } else {
      w_pwin = 0.5 * (1 - cos(2 * M_PI * t / (86400.0 * 365.25)));
      w_psum = 1 - w_pwin;
    }
    w_psum *= (1 - w_midl);
    w_pwin *= (1 - w_midl);
  }

  /* Get maximum SI... */
  s1 = (w_eqn * 63.75 + w_midl * 39.88 + w_psum * 10.73 + w_pwin * 45.58)
    / (w_eqn + w_midl + w_psum + w_pwin);

  /* Scale SI... */
  si *= sia[*n - 1] / s1;

  /* Estimate column density... */
  if (si <= sia[0]) {
    *scd = 0;
    *err = GSL_NAN;
  } else if (si >= sia[*n - 1]) {
    *scd = GSL_POSINF;
    *err = GSL_POSINF;
  } else {
    i = locate_irr(sia, *n, si);
    *scd = LIN(sia[i], scda[i], sia[i + 1], scda[i + 1], si);

    i = locate_irr(sia, *n, si + dsi + 1.0);
    scdp = LIN(sia[i], scda[i], sia[i + 1], scda[i + 1], si + dsi + 1.0);

    i = locate_irr(sia, *n, si - dsi - 1.0);
    scdm = LIN(sia[i], scda[i], sia[i + 1], scda[i + 1], si - dsi - 1.0);

    *err = GSL_MAX(fabs(scdm - *scd), fabs(scdp - *scd));
  }
}
