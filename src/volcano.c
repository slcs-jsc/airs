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
  Detection methods for volcanic ash and sulfur dioxide.
*/

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

  const double ci_nedt = 0.0783, ai_low_bt1_nedt = 0.3698, ai_low_bt2_nedt =
    0.1177, ai_high_bt1_nedt = 0.0766, ai_high_bt2_nedt =
    0.3706, ai_old_bt1_nedt = 0.3155, ai_old_bt2_nedt =
    0.1177, si_high_bt1_nedt = 0.1025, si_high_bt2_nedt =
    0.1373, si_low_bt1_nedt = 0.0799, si_low_bt2_nedt =
    0.0909, si_old_bt1_nedt = 0.1064, si_old_bt2_nedt =
    0.0909, si_oper_bt1_nedt = 0.0884, si_oper_bt2_nedt = 0.1159;

  const int ai_low_nu1 = 641, ai_low_nu2 = 901, ai_high_nu1 =
    1295, ai_high_nu2 = 1162, ai_old_nu1 = 559, ai_old_nu2 = 901, ci_nu =
    1290, si_low_nu1 = 1601, si_low_nu2 = 1526, si_high_nu1 =
    1602, si_high_nu2 = 1551, si_old_nu1 = 1591, si_old_nu2 =
    1526, si_oper_nu1 = 1636, si_oper_nu2 = 1507;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <out.tab> <l1b_file1> [<l1b_file2> ...]");

  /* Create file... */
  printf("Write volcanic emission data: %s\n", argv[1]);
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Loop over HDF files... */
  for (int iarg = 2; iarg < argc; iarg++) {

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
	      "# $9  = ash index (low wavenumbers),"
	      " BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $10 = ash index (low wavenumbers) error [K]\n"
	      "# $11 = ash index (high wavenumbers),"
	      " BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $12 = ash index (high wavenumbers) error [K]\n"
	      "# $13 = ash index (Hoffmann et al., 2014),"
	      " BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $14 = ash index (Hoffmann et al., 2014) error [K]\n",
	      airs_rad_gran.nominal_freq[ci_nu],
	      airs_rad_gran.nominal_freq[ai_low_nu1],
	      airs_rad_gran.nominal_freq[ai_low_nu2],
	      airs_rad_gran.nominal_freq[ai_high_nu1],
	      airs_rad_gran.nominal_freq[ai_high_nu2],
	      airs_rad_gran.nominal_freq[ai_old_nu1],
	      airs_rad_gran.nominal_freq[ai_old_nu2]);
      fprintf(out,
	      "# $15 = SO2 index (low concentrations),"
	      " BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $16 = SO2 index (low concentrations) error [K]\n"
	      "# $17 = SO2 index (high concentrations),"
	      " BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $18 = SO2 index (high concentrations) error [K]\n"
	      "# $19 = SO2 index (operational),"
	      " BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $20 = SO2 index (operational) error [K]\n"
	      "# $21 = SO2 index (Hoffmann et al., 2014),"
	      " BT(%.2f/cm) - BT(%.2f/cm) [K]\n"
	      "# $22 = SO2 index (Hoffmann et al., 2014) error [K]\n",
	      airs_rad_gran.nominal_freq[si_low_nu1],
	      airs_rad_gran.nominal_freq[si_low_nu2],
	      airs_rad_gran.nominal_freq[si_high_nu1],
	      airs_rad_gran.nominal_freq[si_high_nu2],
	      airs_rad_gran.nominal_freq[si_oper_nu1],
	      airs_rad_gran.nominal_freq[si_oper_nu2],
	      airs_rad_gran.nominal_freq[si_old_nu1],
	      airs_rad_gran.nominal_freq[si_old_nu2]);
    }

    /* Flag bad observations... */
    for (int track = 0; track < AIRS_RAD_GEOTRACK; track++)
      for (int xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
	for (int ichan = 0; ichan < AIRS_RAD_CHANNEL; ichan++)
	  if ((airs_rad_gran.state[track][xtrack] != 0)
	      || (airs_rad_gran.ExcludedChans[ichan] > 2)
	      || (airs_rad_gran.CalChanSummary[ichan] & 8)
	      || (airs_rad_gran.CalChanSummary[ichan] & (32 + 64))
	      || (airs_rad_gran.CalFlag[track][ichan] & 16))
	    airs_rad_gran.radiances[track][xtrack][ichan] = GSL_NAN;

    /* Loop over scans... */
    for (int track = 0; track < AIRS_RAD_GEOTRACK; track++) {

      /* Write output... */
      fprintf(out, "\n");

      /* Loop over footprints... */
      for (int xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {

	/* cloud index... */
	double ci =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][ci_nu] * 0.001,
		 airs_rad_gran.nominal_freq[ci_nu]);
	double ci_err =
	  get_noise(ci, ci_nedt, airs_rad_gran.nominal_freq[ci_nu]);

	/* ash index (low wavenumbers)... */
	double ai_low_bt1 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][ai_low_nu1] *
		 0.001, airs_rad_gran.nominal_freq[ai_low_nu1]);
	double ai_low_bt2 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][ai_low_nu2] *
		 0.001, airs_rad_gran.nominal_freq[ai_low_nu2]);
	double ai_low = ai_low_bt1 - ai_low_bt2;
	double ai_low_err =
	  sqrt(gsl_pow_2(get_noise(ai_low_bt1, ai_low_bt1_nedt,
				   airs_rad_gran.nominal_freq[ai_low_nu1]))
	       + gsl_pow_2(get_noise(ai_low_bt2, ai_low_bt2_nedt,
				     airs_rad_gran.nominal_freq
				     [ai_low_nu2])));

	/* ash index (high wavenumbers)... */
	double ai_high_bt1 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][ai_high_nu1] *
		 0.001, airs_rad_gran.nominal_freq[ai_high_nu1]);
	double ai_high_bt2 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][ai_high_nu2] *
		 0.001, airs_rad_gran.nominal_freq[ai_high_nu2]);
	double ai_high = ai_high_bt1 - ai_high_bt2;
	double ai_high_err =
	  sqrt(gsl_pow_2(get_noise(ai_high_bt1, ai_high_bt1_nedt,
				   airs_rad_gran.nominal_freq[ai_high_nu1]))
	       + gsl_pow_2(get_noise(ai_high_bt2, ai_high_bt2_nedt,
				     airs_rad_gran.nominal_freq
				     [ai_high_nu2])));

	/* ash index (old)... */
	double ai_old_bt1 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][ai_old_nu1] *
		 0.001, airs_rad_gran.nominal_freq[ai_old_nu1]);
	double ai_old_bt2 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][ai_old_nu2] *
		 0.001, airs_rad_gran.nominal_freq[ai_old_nu2]);
	double ai_old = ai_old_bt1 - ai_old_bt2;
	double ai_old_err =
	  sqrt(gsl_pow_2(get_noise(ai_old_bt1, ai_old_bt1_nedt,
				   airs_rad_gran.nominal_freq[ai_old_nu1]))
	       + gsl_pow_2(get_noise(ai_old_bt2, ai_old_bt2_nedt,
				     airs_rad_gran.nominal_freq
				     [ai_old_nu2])));

	/* SO2 index (low concentrations)... */
	double si_low_bt1 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][si_low_nu1] *
		 0.001, airs_rad_gran.nominal_freq[si_low_nu1]);
	double si_low_bt2 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][si_low_nu2] *
		 0.001, airs_rad_gran.nominal_freq[si_low_nu2]);
	double si_low = si_low_bt1 - si_low_bt2;
	double si_low_err =
	  sqrt(gsl_pow_2(get_noise(si_low_bt1, si_low_bt1_nedt,
				   airs_rad_gran.nominal_freq[si_low_nu1]))
	       + gsl_pow_2(get_noise(si_low_bt2, si_low_bt2_nedt,
				     airs_rad_gran.nominal_freq
				     [si_low_nu2])));

	/* SO2 index (high concentrations)... */
	double si_high_bt1 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][si_high_nu1] *
		 0.001, airs_rad_gran.nominal_freq[si_high_nu1]);
	double si_high_bt2 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][si_high_nu2] *
		 0.001, airs_rad_gran.nominal_freq[si_high_nu2]);
	double si_high = si_high_bt1 - si_high_bt2;
	double si_high_err =
	  sqrt(gsl_pow_2(get_noise(si_high_bt1, si_high_bt1_nedt,
				   airs_rad_gran.nominal_freq[si_high_nu1]))
	       + gsl_pow_2(get_noise(si_high_bt2, si_high_bt2_nedt,
				     airs_rad_gran.nominal_freq
				     [si_high_nu2])));

	/* SO2 index (operational)... */
	double si_oper_bt1 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][si_oper_nu1] *
		 0.001, airs_rad_gran.nominal_freq[si_oper_nu1]);
	double si_oper_bt2 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][si_oper_nu2] *
		 0.001, airs_rad_gran.nominal_freq[si_oper_nu2]);
	double si_oper = si_oper_bt1 - si_oper_bt2;
	double si_oper_err =
	  sqrt(gsl_pow_2(get_noise(si_oper_bt1, si_oper_bt1_nedt,
				   airs_rad_gran.nominal_freq[si_oper_nu1]))
	       + gsl_pow_2(get_noise(si_oper_bt2, si_oper_bt2_nedt,
				     airs_rad_gran.nominal_freq
				     [si_oper_nu2])));

	/* SO2 index (old)... */
	double si_old_bt1 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][si_old_nu1] *
		 0.001, airs_rad_gran.nominal_freq[si_old_nu1]);
	double si_old_bt2 =
	  BRIGHT(airs_rad_gran.radiances[track][xtrack][si_old_nu2] *
		 0.001, airs_rad_gran.nominal_freq[si_old_nu2]);
	double si_old = si_old_bt1 - si_old_bt2;
	double si_old_err =
	  sqrt(gsl_pow_2(get_noise(si_old_bt1, si_old_bt1_nedt,
				   airs_rad_gran.nominal_freq[si_old_nu1]))
	       + gsl_pow_2(get_noise(si_old_bt2, si_old_bt2_nedt,
				     airs_rad_gran.nominal_freq
				     [si_old_nu2])));

	/* Write output... */
	fprintf(out,
		"%.2f %.4f %.4f %.3f %.4f %.4f %.2f %.2f %.2f %.2f %.2f %.2f "
		"%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",
		airs_rad_gran.Time[track][xtrack] - 220838400,
		airs_rad_gran.Longitude[track][xtrack],
		airs_rad_gran.Latitude[track][xtrack],
		airs_rad_gran.satheight[track],
		airs_rad_gran.sat_lon[track],
		airs_rad_gran.sat_lat[track],
		ci, ci_err, ai_low, ai_low_err, ai_high, ai_high_err, ai_old,
		ai_old_err, si_low, si_low_err, si_high, si_high_err, si_oper,
		si_oper_err, si_old, si_old_err);
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

  const double nesr = NESR(250.0, dt250, nu);
  
  return NEDT(bt, nesr, nu);
}
