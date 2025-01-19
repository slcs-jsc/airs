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
  Extract maps of gravity wave perturbations.
*/

#include "libairs.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Fill data gaps in perturbation data. */
double fill_array(
  double var[PERT_NTRACK][PERT_NXTRACK],
  int ntrack,
  int itrack,
  int ixtrack);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert, *pert2;
  static wave_t wave;

  char set[LEN], pertname[LEN];

  double nedt = 0, sza2 = 0;

  int orb = 0;

  FILE *out;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <pert.nc> <map.tab>");

  /* Get control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  const int bg_poly_x =
    (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "0", NULL);
  const int bg_poly_y =
    (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  const int bg_smooth_x =
    (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  const int bg_smooth_y =
    (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "0", NULL);
  const double gauss_fwhm = scan_ctl(argc, argv, "GAUSS_FWHM", -1, "0", NULL);
  const int ham_iter = (int) scan_ctl(argc, argv, "HAM_ITER", -1, "0", NULL);
  const int med_dx = (int) scan_ctl(argc, argv, "MED_DX", -1, "0", NULL);
  const double var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "0", NULL);
  scan_ctl(argc, argv, "SET", -1, "full", set);
  const int orbit = (int) scan_ctl(argc, argv, "ORBIT", -1, "-999", NULL);
  const double orblat = scan_ctl(argc, argv, "ORBLAT", -1, "0", NULL);
  const double t0 = scan_ctl(argc, argv, "T0", -1, "-1e100", NULL);
  const double t1 = scan_ctl(argc, argv, "T1", -1, "1e100", NULL);
  const double sza0 = scan_ctl(argc, argv, "SZA0", -1, "-1e100", NULL);
  const double sza1 = scan_ctl(argc, argv, "SZA1", -1, "1e100", NULL);
  const double dt230 = scan_ctl(argc, argv, "DT230", -1, "0.16", NULL);
  const double nu = scan_ctl(argc, argv, "NU", -1, "2345.0", NULL);
  const int fill = (int) scan_ctl(argc, argv, "FILL", -1, "0", NULL);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);
  ALLOC(pert2, pert_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], pertname, pert);

  /* Recalculate background and perturbations... */
  if (bg_poly_x > 0 || bg_poly_y > 0 ||
      bg_smooth_x > 0 || bg_smooth_y > 0 ||
      gauss_fwhm > 0 || ham_iter > 0 || med_dx > 0 || var_dh > 0) {

    /* Convert to wave analysis struct... */
    pert2wave(pert, &wave, 0, pert->ntrack - 1, 0, pert->nxtrack - 1);

    /* Estimate background... */
    background_poly(&wave, bg_poly_x, bg_poly_y);
    background_smooth(&wave, bg_smooth_x, bg_smooth_y);

    /* Gaussian filter... */
    gauss(&wave, gauss_fwhm);

    /* Hamming filter... */
    hamming(&wave, ham_iter);

    /* Median filter... */
    median(&wave, med_dx);

    /* Compute variance... */
    variance(&wave, var_dh);

    /* Copy data... */
    for (int ix = 0; ix < wave.nx; ix++)
      for (int iy = 0; iy < wave.ny; iy++) {
	pert->pt[iy][ix] = wave.pt[ix][iy];
	pert->var[iy][ix] = wave.var[ix][iy];
      }
  }

  /* Fill data gaps... */
  if (fill)
    for (int itrack = 0; itrack < pert->ntrack; itrack++)
      for (int ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {
	if (!gsl_finite(pert->dc[itrack][ixtrack]))
	  pert->dc[itrack][ixtrack]
	    = fill_array(pert->dc, pert->ntrack, itrack, ixtrack);
	if (!gsl_finite(pert->bt[itrack][ixtrack]))
	  pert->bt[itrack][ixtrack]
	    = fill_array(pert->bt, pert->ntrack, itrack, ixtrack);
	if (!gsl_finite(pert->pt[itrack][ixtrack]))
	  pert->pt[itrack][ixtrack]
	    = fill_array(pert->pt, pert->ntrack, itrack, ixtrack);
	if (!gsl_finite(pert->var[itrack][ixtrack]))
	  pert->var[itrack][ixtrack]
	    = fill_array(pert->var, pert->ntrack, itrack, ixtrack);
      }

  /* Interpolate to fine grid... */
  memcpy(pert2, pert, sizeof(pert_t));

  /* Create output file... */
  printf("Write perturbation data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = solar zenith angle [deg]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = 8mu brightness temperature [K]\n"
	  "# $6 = %s brightness temperature [K]\n"
	  "# $7 = %s brightness temperature perturbation [K]\n"
	  "# $8 = %s brightness temperature variance [K^2]\n"
	  "# $9 = along-track index\n"
	  "# $10 = across-track index\n", pertname, pertname, pertname);

  /* Write data... */
  for (int itrack = 0; itrack < pert->ntrack; itrack++) {

    /* Count orbits... */
    if (itrack > 0)
      if (pert->lat[itrack - 1][pert->nxtrack / 2] <= orblat
	  && pert->lat[itrack][pert->nxtrack / 2] >= orblat)
	orb++;

    /* Write output... */
    fprintf(out, "\n");

    /* Check for data gaps... */
    if (itrack > 0 && pert->time[itrack][pert->nxtrack / 2]
	- pert->time[itrack - 1][pert->nxtrack / 2] >= 10)
      fprintf(out, "\n");

    /* Loop over scan... */
    for (int ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {

      /* Check data... */
      if (pert->lon[itrack][ixtrack] < -180
	  || pert->lon[itrack][ixtrack] > 180
	  || pert->lat[itrack][ixtrack] < -90
	  || pert->lat[itrack][ixtrack] > 90)
	continue;

      /* Get ascending/descending flag... */
      int asc =
	(pert->lat[itrack > 0 ? itrack : itrack + 1][pert->nxtrack / 2]
	 > pert->lat[itrack > 0 ? itrack - 1 : itrack][pert->nxtrack / 2]);

      /* Calculate solar zenith angle... */
      if (sza0 >= -1e10 && sza0 <= 1e10 && sza1 >= -1e10 && sza1 <= 1e10)
	sza2 = sza(pert->time[itrack][ixtrack], pert->lon[itrack][ixtrack],
		   pert->lat[itrack][ixtrack]);

      /* Estimate noise... */
      if (dt230 > 0) {
	const double nesr = PLANCK(230.0 + dt230, nu) - PLANCK(230.0, nu);
	const double tbg =
	  pert->bt[itrack][ixtrack] - pert->pt[itrack][ixtrack];
	nedt = BRIGHT(PLANCK(tbg, nu) + nesr, nu) - tbg;
      }

      /* Write data... */
      if (orbit < 0 || orb == orbit)
	if (set[0] == 'f' || (set[0] == 'a' && asc)
	    || (set[0] == 'd' && !asc))
	  if (pert->time[itrack][ixtrack] >= t0
	      && pert->time[itrack][ixtrack] <= t1
	      && sza2 >= sza0 && sza2 <= sza1)
	    fprintf(out, "%.2f %g %g %g %g %g %g %g %d %d\n",
		    pert->time[itrack][ixtrack],
		    sza(pert->time[itrack][ixtrack],
			pert->lon[itrack][ixtrack],
			pert->lat[itrack][ixtrack]),
		    pert->lon[itrack][ixtrack], pert->lat[itrack][ixtrack],
		    pert->dc[itrack][ixtrack], pert->bt[itrack][ixtrack],
		    pert->pt[itrack][ixtrack],
		    pert->var[itrack][ixtrack] - gsl_pow_2(nedt), itrack,
		    ixtrack);
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(pert);
  free(pert2);

  return EXIT_SUCCESS;
}

/************************************************************************/

double fill_array(
  double var[PERT_NTRACK][PERT_NXTRACK],
  int ntrack,
  int itrack,
  int ixtrack) {

  double d1 = 0, d2 = 0, v1 = 0, v2 = 0;

  /* Find nearest neighbours... */
  for (int i = itrack + 1; i < ntrack; i++)
    if (gsl_finite(var[i][ixtrack])) {
      d1 = fabs((double) i - (double) itrack);
      v1 = var[i][ixtrack];
      break;
    }
  for (int i = itrack - 1; i >= 0; i--)
    if (gsl_finite(var[i][ixtrack])) {
      d2 = fabs((double) i - (double) itrack);
      v2 = var[i][ixtrack];
      break;
    }

  /* Interpolate... */
  if (d1 + d2 > 0)
    return (d2 * v1 + d1 * v2) / (d1 + d2);
  else
    return GSL_NAN;
}
