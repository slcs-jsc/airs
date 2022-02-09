#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static wave_t *wave;

  FILE *out;

  char method[LEN], filename[LEN];

  double Amax, phimax, lhmax, kxmax, kymax, alphamax, betamax, chisq, lh,
    alpha;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <spec.tab>");

  /* Get control parameters... */
  int nx = (int) scan_ctl(argc, argv, "NX", -1, "256", NULL);
  int ny = (int) scan_ctl(argc, argv, "NY", -1, "256", NULL);
  double sx = scan_ctl(argc, argv, "SX", -1, "1000.0", NULL);
  double sy = scan_ctl(argc, argv, "SY", -1, "1000.0", NULL);
  double amp = scan_ctl(argc, argv, "AMP", -1, "1.0", NULL);
  double phi = scan_ctl(argc, argv, "PHI", -1, "0.0", NULL);
  double lx0 = scan_ctl(argc, argv, "LX0", -1, "100.0", NULL);
  double lx1 = scan_ctl(argc, argv, "LX1", -1, "100.0", NULL);
  double dlx = scan_ctl(argc, argv, "DLX", -1, "10.0", NULL);
  double ly0 = scan_ctl(argc, argv, "LY0", -1, "200.0", NULL);
  double ly1 = scan_ctl(argc, argv, "LY1", -1, "200.0", NULL);
  double dly = scan_ctl(argc, argv, "DLY", -1, "10.0", NULL);
  double fwhm = scan_ctl(argc, argv, "FWHM", -1, "0.0", NULL);
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
  ALLOC(wave, wave_t, 1);

  /* Create output file... */
  printf("Write spectral data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = real amplitude [K]\n"
	  "# $2 = real phase [deg]\n"
	  "# $3 = real horizontal wavelength [km]\n"
	  "# $4 = real wavelength x [km]\n"
	  "# $5 = real wavelength y [km]\n"
	  "# $6 = real rotation angle alpha [deg]\n");
  fprintf(out,
	  "# $7 = amplitude [K]\n"
	  "# $8 = phase [deg]\n"
	  "# $9 = horizontal wavelength [km]\n"
	  "# $10 = wavelength x [km]\n"
	  "# $11 = wavelength y [km]\n"
	  "# $12 = rotation angle alpha [deg]\n"
	  "# $13 = rotation angle beta [deg]\n"
	  "# $14 = chi^2 of fit [K^2]\n"
	  "# $15 = box size x\n" "# $16 = box size y\n");

  /* Set grid... */
  wave->nx = nx;
  wave->ny = ny;
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {
      wave->x[ix] = sx / (nx - 1.) * ix - sx / 2.;
      wave->y[iy] = sy / (ny - 1.) * iy - sy / 2.;
    }

  /* Loop over wavelengths... */
  for (double lx = lx0; lx <= lx1; lx += dlx) {

    /* Write output... */
    fprintf(out, "\n");

    /* Loop over wavelengths... */
    for (double ly = ly0; ly <= ly1; ly += dly) {

      /* Write info... */
      printf("Analyze lx = %g km and ly = %g km...\n", lx, ly);

      /* Get horizontal wavelength... */
      lh =
	2 * M_PI / sqrt(gsl_pow_2(2 * M_PI / lx) + gsl_pow_2(2 * M_PI / ly));

      /* Get propagation direction in xy-plane... */
      alpha = 90. - 180. / M_PI * atan2(2 * M_PI / lx, 2 * M_PI / ly);

      /* Create wave field... */
      create_wave(wave, amp, lx, ly, phi, fwhm);

      /* Estimate background... */
      background_poly(wave, bg_poly_x, bg_poly_y);
      background_smooth(wave, bg_smooth_x, bg_smooth_y);

      /* Compute variance... */
      variance(wave, var_dh);

      /* Interpolate to regular grid... */
      intpol_x(wave, inter_x);

      /* Get wave characteristics... */
      if (method[0] == 'p' || method[0] == 'P') {
	sprintf(filename, "period_%g_%g.tab", lx, ly);
	period(wave, lxymax, dlxy, &Amax, &phimax, &lhmax, &kxmax, &kymax,
	       &alphamax, &betamax, output ? filename : NULL);
      } else if (method[0] == 'f' || method[0] == 'F') {
	sprintf(filename, "fft_%g_%g.tab", lx, ly);
	fft(wave, &Amax, &phimax, &lhmax, &kxmax, &kymax,
	    &alphamax, &betamax, output ? filename : NULL);
      } else
	ERRMSG("Unkown method!");

      /* Evaluate fit... */
      fit_wave(wave, Amax, phimax, kxmax, kymax, &chisq);

      /* Save wave struct... */
      sprintf(filename, "wave_%g_%g.tab", lx, ly);
      if (output)
	write_wave(filename, wave);

      /* Write results... */
      fprintf(out, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %d\n",
	      amp, phi, lh, lx, ly, alpha, Amax, phimax, lhmax,
	      2 * M_PI / kxmax, 2 * M_PI / kymax, alphamax, betamax,
	      chisq, wave->nx, wave->ny);
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(wave);

  return EXIT_SUCCESS;
}
