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
  Spectral analysis of gravity wave perturbations.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static wave_t *wave;
  static pert_t *pert;

  FILE *out;

  char method[LEN], filename[LEN], pertname[LEN];

  double Amax, phimax, lhmax, kxmax, kymax, alphamax, betamax, chisq;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <pert.nc> <spec.tab>");

  /* Get control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  int track0 = (int) scan_ctl(argc, argv, "TRACK0", -1, "", NULL);
  int track1 = (int) scan_ctl(argc, argv, "TRACK1", -1, "", NULL);
  int xtrack0 = (int) scan_ctl(argc, argv, "XTRACK0", -1, "0", NULL);
  int xtrack1 = (int) scan_ctl(argc, argv, "XTRACK1", -1, "89", NULL);
  int strack = (int) scan_ctl(argc, argv, "STRACK", -1, "15", NULL);
  int sxtrack = (int) scan_ctl(argc, argv, "SXTRACK", -1, "15", NULL);
  int dtrack = (int) scan_ctl(argc, argv, "DTRACK", -1, "5", NULL);
  int dxtrack = (int) scan_ctl(argc, argv, "DXTRACK", -1, "5", NULL);
  int inter_x = (int) scan_ctl(argc, argv, "INTER_X", -1, "0", NULL);
  int bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "5", NULL);
  int bg_poly_y = (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  int bg_smooth_x = (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  int bg_smooth_y = (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "0", NULL);
  double var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "100.0", NULL);
  double lxymax = scan_ctl(argc, argv, "LXYMAX", -1, "1000.0", NULL);
  double dlxy = scan_ctl(argc, argv, "DLXY", -1, "10.0", NULL);
  scan_ctl(argc, argv, "METHOD", -1, "P", method);
  int output = (int) scan_ctl(argc, argv, "OUTPUT", -1, "1", NULL);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);
  ALLOC(wave, wave_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], pertname, pert);

  /* Convert to wave analysis struct... */
  pert2wave(pert, wave, 0, pert->ntrack - 1, 0, pert->nxtrack - 1);

  /* Estimate background... */
  background_poly(wave, bg_poly_x, bg_poly_y);
  background_smooth(wave, bg_smooth_x, bg_smooth_y);

  /* Compute variance... */
  variance(wave, var_dh);

  /* Copy data... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {
      pert->pt[iy][ix] = wave->pt[ix][iy];
      pert->var[iy][ix] = wave->var[ix][iy];
    }

  /* Interpolate to regular grid... */
  intpol_x(wave, inter_x);

  /* Create output file... */
  printf("Write spectral data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = longitude [deg]\n"
	  "# $3 = latitude [deg]\n"
	  "# $4 = amplitude [K]\n"
	  "# $5 = phase [deg]\n"
	  "# $6 = horizontal wavelength [km]\n"
	  "# $7 = wavenumber x [1/km]\n"
	  "# $8 = wavenumber y [1/km]\n"
	  "# $9 = rotation angle alpha [deg]\n"
	  "# $10 = rotation angle beta [deg]\n");
  fprintf(out,
	  "# $11 = chi^2 of fit [K^2]\n"
	  "# $12 = along-track index minimum\n"
	  "# $13 = along-track index maximum\n"
	  "# $14 = across-track index minimum\n"
	  "# $15 = across-track index maximum\n"
	  "# $16 = box size x\n" "# $17 = box size y\n\n");

  /* Loop over blocks... */
  for (int track = track0; track <= track1 - strack + 1; track += dtrack)
    for (int xtrack = xtrack0; xtrack <= xtrack1 - sxtrack + 1;
	 xtrack += dxtrack) {

      /* Write info... */
      printf("Analyze track = %d and xtrack = %d ...\n", track, xtrack);

      /* Convert to wave analysis struct... */
      pert2wave(pert, wave, track, track + strack - 1,
		xtrack, xtrack + sxtrack - 1);

      /* Get wave characteristics... */
      if (method[0] == 'p' || method[0] == 'P') {
	sprintf(filename, "period_%d_%d.tab", track, xtrack);
	period(wave, lxymax, dlxy, &Amax, &phimax, &lhmax, &kxmax, &kymax,
	       &alphamax, &betamax, output ? filename : NULL);
      } else if (method[0] == 'f' || method[0] == 'F') {
	sprintf(filename, "fft_%d_%d.tab", track, xtrack);
	fft(wave, &Amax, &phimax, &lhmax, &kxmax, &kymax,
	    &alphamax, &betamax, output ? filename : NULL);
      } else
	ERRMSG("Unkown method!");

      /* Evaluate fit... */
      fit_wave(wave, Amax, phimax, kxmax, kymax, &chisq);

      /* Save wave struct... */
      if (output) {
	sprintf(filename, "wave_%d_%d.tab", track, xtrack);
	write_wave(filename, wave);
      }

      /* Write results... */
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %d %d %d %d %d %d\n",
	      pert->time[track + strack / 2][xtrack + sxtrack / 2],
	      pert->lon[track + strack / 2][xtrack + sxtrack / 2],
	      pert->lat[track + strack / 2][xtrack + sxtrack / 2],
	      Amax, phimax, lhmax, kxmax, kymax, alphamax, betamax,
	      chisq, track, track + strack - 1,
	      xtrack, xtrack + sxtrack - 1, wave->nx, wave->ny);
    }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(pert);
  free(wave);

  return EXIT_SUCCESS;
}
