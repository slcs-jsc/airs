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
  Identify gravity wave events in AIRS data.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  static wave_t *wave;

  static FILE *out;

  static char pertname[LEN];

  const int dtrack = 15, dxtrack = 15;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <events.tab> <pert1.nc> [<pert2.nc> ...]");

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
  const double var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "0", NULL);
  const double varmin = scan_ctl(argc, argv, "VARMIN", -1, "", NULL);
  const double dt230 = scan_ctl(argc, argv, "DT230", -1, "0.16", NULL);
  const double nu = scan_ctl(argc, argv, "NU", -1, "2345.0", NULL);

  /* Alloc... */
  ALLOC(pert, pert_t, 1);

  /* Create file... */
  printf("Write event data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = longitude [deg]\n"
	  "# $3 = latitude [deg]\n" "# $4 = maximum variance [K^2]\n\n");

  /* Loop over perturbation files... */
  for (int iarg = 3; iarg < argc; iarg++) {

    /* Read perturbation data... */
    FILE *in;
    if (!(in = fopen(argv[iarg], "r")))
      continue;
    else {
      fclose(in);
      read_pert(argv[iarg], pertname, pert);
    }

    /* Recalculate background and perturbations... */
    if (bg_poly_x > 0 || bg_poly_y > 0 ||
	bg_smooth_x > 0 || bg_smooth_y > 0 || gauss_fwhm > 0 || var_dh > 0) {

      /* Allocate... */
      ALLOC(wave, wave_t, 1);

      /* Convert to wave analysis struct... */
      pert2wave(pert, wave, 0, pert->ntrack - 1, 0, pert->nxtrack - 1);

      /* Estimate background... */
      background_poly(wave, bg_poly_x, bg_poly_y);
      background_smooth(wave, bg_smooth_x, bg_smooth_y);

      /* Gaussian filter... */
      gauss(wave, gauss_fwhm);

      /* Compute variance... */
      variance(wave, var_dh);

      /* Copy data... */
      for (int ix = 0; ix < wave->nx; ix++)
	for (int iy = 0; iy < wave->ny; iy++) {
	  pert->pt[iy][ix] = wave->pt[ix][iy];
	  pert->var[iy][ix] = wave->var[ix][iy];
	}

      /* Free... */
      free(wave);
    }

    /* Apply noise correction... */
    if (dt230 > 0)
      for (int itrack = 0; itrack < pert->ntrack; itrack++)
	for (int ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {
	  const double nesr = PLANCK(230.0 + dt230, nu) - PLANCK(230.0, nu);
	  const double tbg =
	    pert->bt[itrack][ixtrack] - pert->pt[itrack][ixtrack];
	  const double nedt = BRIGHT(PLANCK(tbg, nu) + nesr, nu) - tbg;
	  pert->var[itrack][ixtrack] -= gsl_pow_2(nedt);
	}

    /* Find local maxima... */
    for (int itrack = 0; itrack < pert->ntrack; itrack += 2 * dtrack)
      for (int ixtrack = dxtrack / 2; ixtrack < pert->nxtrack;
	   ixtrack += 2 * dxtrack) {

	/* Init... */
	double varmax = 0;
	int itrackmax = -999;
	int ixtrackmax = -999;

	/* Loop over box... */
	for (int itrack2 = itrack;
	     itrack2 < GSL_MIN(itrack + dtrack, pert->ntrack); itrack2++)
	  for (int ixtrack2 = ixtrack;
	       ixtrack2 < GSL_MIN(ixtrack + dxtrack, pert->nxtrack);
	       ixtrack2++)
	    if (pert->var[itrack2][ixtrack2] >= varmax) {
	      varmax = pert->var[itrack2][ixtrack2];
	      itrackmax = itrack2;
	      ixtrackmax = ixtrack2;
	    }

	/* Report event... */
	if (itrackmax >= 0 && ixtrackmax >= 0 && varmax >= varmin)
	  fprintf(out, "%.2f %g %g %g\n",
		  pert->time[itrackmax][ixtrackmax],
		  pert->lon[itrackmax][ixtrackmax],
		  pert->lat[itrackmax][ixtrackmax],
		  pert->var[itrackmax][ixtrackmax]);
      }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}
