#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static ret_t ret;
  static wave_t wave, wave2;

  FILE *out;

  double mu, mu2, nedt, nedt2;

  int ip;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <airs.nc> <noise.tab>");

  /* Read AIRS data... */
  read_retr(argv[2], &ret);

  /* Create file... */
  printf("Write noise data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = altitude [km]\n"
	  "# $2 = longitude [deg]\n"
	  "# $3 = latitude [deg]\n"
	  "# $4 = mean temperature (retrieval) [K]\n"
	  "# $5 = noise estimate (retrieval) [K]\n"
	  "# $6 = mean temperature (a priori) [K]\n"
	  "# $7 = noise estimate (a priori) [K]\n\n");

  /* Loop over altitudes... */
  for (ip = 0; ip < ret.np; ip++) {

    /* Convert retrieval data to wave struct... */
    ret2wave(&ret, &wave, 1, ip);
    ret2wave(&ret, &wave2, 2, ip);

    /* Estimate noise... */
    noise(&wave, &mu, &nedt);
    noise(&wave2, &mu2, &nedt2);

    /* Estimate noise... */
    fprintf(out, "%g %g %g %g %g %g %g\n",
	    wave.z,
	    wave.lon[wave.nx / 2][wave.ny / 2],
	    wave.lat[wave.nx / 2][wave.ny / 2], mu, nedt, mu2, nedt2);
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
