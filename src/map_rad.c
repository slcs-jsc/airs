#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;
  static wave_t wave, wave2;

  double gauss_fwhm, nu, var_dh;

  int bg_poly_x, bg_poly_y, bg_smooth_x, bg_smooth_y;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <l1b_file1> <l1b_file2> <nu> <wave.tab>");

  /* Get control parameters... */
  bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "5", NULL);
  bg_poly_y = (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  bg_smooth_x = (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  bg_smooth_y = (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "0", NULL);
  gauss_fwhm = scan_ctl(argc, argv, "GAUSS_FWHM", -1, "0", NULL);
  var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "0", NULL);

  /* Get channel.. */
  nu = atof(argv[4]);

  /* Read AIRS data... */
  printf("Read AIRS Level-1B data file: %s\n", argv[2]);
  airs_rad_rdr(argv[2], &airs_rad_gran);

  /* Convert radiance data to wave struct... */
  rad2wave(&airs_rad_gran, &nu, 1, &wave);

  /* Check if second file is available... */
  if (argv[3][0] != '-') {

    /* Read AIRS data... */
    printf("Read AIRS Level-1B data file: %s\n", argv[3]);
    airs_rad_rdr(argv[3], &airs_rad_gran);

    /* Convert radiance data to wave struct... */
    rad2wave(&airs_rad_gran, &nu, 1, &wave2);

    /* Merge with first file... */
    merge_y(&wave, &wave2);
  }

  /* Compute background... */
  background_poly(&wave, bg_poly_x, bg_poly_y);
  background_smooth(&wave, bg_smooth_x, bg_smooth_y);

  /* Gaussian filter... */
  gauss(&wave, gauss_fwhm);

  /* Compute variance... */
  variance(&wave, var_dh);

  /* Write files... */
  write_wave(argv[5], &wave);

  return EXIT_SUCCESS;
}
