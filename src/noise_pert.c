#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;
  static wave_t wave;

  FILE *out;

  char pertname[LEN];

  double maxvar, mu, nedt = -1e99, nedt_old;

  int bsize, itrack;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <pert.nc> <noise.tab>");

  /* Read control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  bsize = (int) scan_ctl(argc, argv, "BSIZE", -1, "-999", NULL);
  maxvar = (int) scan_ctl(argc, argv, "MAXVAR", -1, "-999", NULL);

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
  for (itrack = 0; itrack < pert->ntrack; itrack += bsize) {

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
