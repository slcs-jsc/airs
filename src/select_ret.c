#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static ret_t ret;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> [<airs1.nc> <airs2.c> ...]");

  /* Get control parameters... */
  double lon0 = (int) scan_ctl(argc, argv, "SELECT_LON0", -1, "-180", NULL);
  double lon1 = (int) scan_ctl(argc, argv, "SELECT_LON1", -1, "180", NULL);
  double lat0 = (int) scan_ctl(argc, argv, "SELECT_LAT0", -1, "-90", NULL);
  double lat1 = (int) scan_ctl(argc, argv, "SELECT_LAT1", -1, "90", NULL);

  /* Loop over retrieval files... */
  for (int i = 2; i < argc; i++) {

    /* Read AIRS data... */
    read_retr(argv[i], &ret);

    /* Check position... */
    if (ret.lon[ret.nds / 2][0] >= lon0
	&& ret.lon[ret.nds / 2][0] <= lon1
	&& ret.lat[ret.nds / 2][0] >= lat0 && ret.lat[ret.nds / 2][0] <= lat1)
      printf("select: %s %.2f %g %g\n", argv[i], ret.time[ret.nds / 2][0],
	     ret.lon[ret.nds / 2][0], ret.lat[ret.nds / 2][0]);
  }

  return EXIT_SUCCESS;
}
