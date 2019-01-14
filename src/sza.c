#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  double jsec, lon, lat;

  /* Check arguments... */
  if (argc != 4)
    ERRMSG("Give parameters: <jsec> <lon> <lat>");

  /* Read arguments... */
  jsec = atof(argv[1]);
  lon = atof(argv[2]);
  lat = atof(argv[3]);

  /* Compute solar zenith angle... */
  printf("%g\n", sza(jsec, lon, lat));

  return EXIT_SUCCESS;
}
