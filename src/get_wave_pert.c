#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static wave_t wave;
  static pert_t *pert;

  char method[LEN], pertname[LEN];

  double var_dh, Amax, phimax, lhmax, alphamax, betamax;

  int bg_poly_x, bg_poly_y, bg_smooth_x, bg_smooth_y, inter_x,
    dtrack, dxtrack, track0, xtrack0;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <pert.nc>");

  /* Get control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  track0 = (int) scan_ctl(argc, argv, "TRACK0", -1, "", NULL);
  xtrack0 = (int) scan_ctl(argc, argv, "XTRACK0", -1, "", NULL);
  dtrack = (int) scan_ctl(argc, argv, "DTRACK", -1, "20", NULL);
  dxtrack = (int) scan_ctl(argc, argv, "DXTRACK", -1, "20", NULL);
  inter_x = (int) scan_ctl(argc, argv, "INTER_X", -1, "0", NULL);
  bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "5", NULL);
  bg_poly_y = (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  bg_smooth_x = (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  bg_smooth_y = (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "7", NULL);
  var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "100", NULL);
  scan_ctl(argc, argv, "METHOD", -1, "P", method);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], pertname, pert);

  /* Check indices... */
  if (track0 < 0 || track0 >= pert->ntrack)
    ERRMSG("Along-track index out of range!");
  if (xtrack0 < 0 || xtrack0 >= pert->nxtrack)
    ERRMSG("Across-track index out of range!");

  /* Convert to wave analysis struct... */
  pert2wave(pert, &wave,
	    track0 - dtrack, track0 + dtrack,
	    xtrack0 - dxtrack, xtrack0 + dxtrack);

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

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}
