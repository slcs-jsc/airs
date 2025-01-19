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
  Extract maps of radiance data.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;
  static wave_t wave, wave2;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <l1b_file1> <l1b_file2> <nu> <wave.tab>");

  /* Get control parameters... */
  const int bg_poly_x =
    (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "5", NULL);
  const int bg_poly_y =
    (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  const int bg_smooth_x =
    (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  const int bg_smooth_y =
    (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "0", NULL);
  const double gauss_fwhm = scan_ctl(argc, argv, "GAUSS_FWHM", -1, "0", NULL);
  const double var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "0", NULL);

  /* Get channel.. */
  double nu = atof(argv[4]);

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
