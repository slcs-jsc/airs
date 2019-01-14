#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static wave_t wave;

  char method[LEN];

  double amp, dx, dy, lx, ly, phi, fwhm, var_dh,
    nedt, Amax, phimax, lhmax, alphamax, betamax;

  int bg_poly_x, bg_poly_y, bg_smooth_x, bg_smooth_y, inter_x, ix, iy, nx, ny;

  /* Check arguments... */
  if (argc < 2)
    ERRMSG("Give parameters: <ctl>");

  /* Get control parameters... */
  nx = (int) scan_ctl(argc, argv, "NX", -1, "90", NULL);
  ny = (int) scan_ctl(argc, argv, "NY", -1, "135", NULL);
  dx = scan_ctl(argc, argv, "DX", -1, "18", NULL);
  dy = scan_ctl(argc, argv, "DY", -1, "18", NULL);
  amp = scan_ctl(argc, argv, "AMP", -1, "1", NULL);
  phi = scan_ctl(argc, argv, "PHI", -1, "0", NULL);
  lx = scan_ctl(argc, argv, "LX", -1, "0", NULL);
  ly = scan_ctl(argc, argv, "LY", -1, "0", NULL);
  fwhm = scan_ctl(argc, argv, "FWHM", -1, "0", NULL);
  nedt = scan_ctl(argc, argv, "NOISE", -1, "0", NULL);
  inter_x = (int) scan_ctl(argc, argv, "INTER_X", -1, "0", NULL);
  bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "5", NULL);
  bg_poly_y = (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  bg_smooth_x = (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  bg_smooth_y = (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "7", NULL);
  var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "100", NULL);
  scan_ctl(argc, argv, "METHOD", -1, "P", method);

  /* Set grid... */
  wave.nx = nx;
  wave.ny = ny;
  for (ix = 0; ix < nx; ix++)
    wave.x[ix] = (ix - nx / 2) * dx;
  for (iy = 0; iy < ny; iy++)
    wave.y[iy] = (iy - ny / 2) * dy;

  /* Init wave... */
  create_background(&wave);
  create_wave(&wave, amp, lx, ly, phi, fwhm);
  create_noise(&wave, nedt);

  /* Interpolate to regular grid... */
  intpol_x(&wave, inter_x);

  /* Estimate background... */
  background_poly(&wave, bg_poly_x, bg_poly_y);
  background_smooth(&wave, bg_smooth_x, bg_smooth_y);

  /* Compute variance... */
  variance(&wave, var_dh);

  /* Get wave characteristics... */
  if (method[0] == 'p' || method[0] == 'P')
    period(&wave, &Amax, &phimax, &lhmax, &alphamax, &betamax, "period.tab");
  if (method[0] == 'f' || method[0] == 'F')
    fft(&wave, &Amax, &phimax, &lhmax, &alphamax, &betamax, "period.tab");

  /* Save wave struct... */
  write_wave("wave.tab", &wave);

  /* Write results... */
  PRINT("%g", Amax);
  PRINT("%g", phimax);
  PRINT("%g", lhmax);
  PRINT("%g", alphamax);
  PRINT("%g", betamax);

  return EXIT_SUCCESS;
}
