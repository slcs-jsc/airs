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
  Estimate noise based on perturbation data.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;
  static wave_t wave;

  FILE *out;

  char pertname[LEN];

  double mu, nedt = -1e99, nedt_old;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <pert.nc> <noise.tab>");

  /* Read control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  int bsize = (int) scan_ctl(argc, argv, "BSIZE", -1, "-999", NULL);
  const double maxvar = scan_ctl(argc, argv, "MAXVAR", -1, "-999", NULL);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], pertname, pert);

  /* Set block size... */
  if (bsize < 0)
    bsize = pert->nxtrack;

  /* Create file... */
  printf("Write noise data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = longitude [deg]\n"
	  "# $2 = latitude [deg]\n"
	  "# $3 = mean brightness temperature [K]\n"
	  "# $4 = noise estimate [K]\n\n");

  /* Loop over granules... */
  for (int itrack = 0; itrack < pert->ntrack; itrack += bsize) {

    /* Convert retrieval data to wave struct... */
    pert2wave(pert, &wave, itrack, itrack + bsize,
	      pert->nxtrack / 2 - bsize / 2, pert->nxtrack / 2 + bsize / 2);

    /* Estimate noise... */
    nedt_old = nedt;
    noise(&wave, &mu, &nedt);

    /* Write output... */
    if (maxvar <= 0
	|| fabs(200 * (nedt - nedt_old) / (nedt + nedt_old)) < maxvar)
      fprintf(out, "%g %g %g %g\n", wave.lon[wave.nx / 2][wave.ny / 2],
	      wave.lat[wave.nx / 2][wave.ny / 2], mu, nedt);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}
