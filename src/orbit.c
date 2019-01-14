#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;

  FILE *out;

  int i, track, xtrack;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG
      ("Give parameters: <orbit.tab> <airs_l1b_file> [ <airs_l1b_file2> ... ]");

  /* Create file... */
  printf("Write orbit data: %s\n", argv[1]);
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = satellite longitude [deg]\n"
	  "# $3 = satellite latitude [deg]\n"
	  "# $4 = footprint longitude [deg]\n"
	  "# $5 = footprint latitude [deg]\n");

  /* Loop over files... */
  for (i = 2; i < argc; i++) {

    /* Read AIRS data... */
    printf("Read AIRS Level-1B data file: %s\n", argv[i]);
    airs_rad_rdr(argv[i], &airs_rad_gran);

    /* Write data... */
    for (track = 0; track < AIRS_RAD_GEOTRACK; track++) {
      fprintf(out, "\n");
      for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
	fprintf(out, "%.2f %g %g %g %g\n",
		airs_rad_gran.Time[track][xtrack] - 220838400,
		airs_rad_gran.sat_lon[track],
		airs_rad_gran.sat_lat[track],
		airs_rad_gran.Longitude[track][xtrack],
		airs_rad_gran.Latitude[track][xtrack]);
    }
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
