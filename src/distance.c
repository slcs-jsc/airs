#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  double lat0, lat1, lon0, lon1, x0[3], x1[3];

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <lon0> <lat0> <lon1> <lat1>");

  /* Read geolocations... */
  lon0 = atof(argv[1]);
  lat0 = atof(argv[2]);
  lon1 = atof(argv[3]);
  lat1 = atof(argv[4]);

  /* Write distance to stdout... */
  geo2cart(0, lon0, lat0, x0);
  geo2cart(0, lon1, lat1, x1);
  printf("%g\n", DIST(x0, x1));

  return EXIT_SUCCESS;
}
