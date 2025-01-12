#include "libairs.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Write results to file. */
void write_results(
  FILE * out,
  pert_t * pert,
  int track0,
  int xtrack0,
  int orb,
  double dmin,
  double obsz);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  FILE *out;

  char pertname[LEN];

  double dmin = 1e100, lon0, lat0, orblat, rmax, obsz, x0[3], x1[3];

  int orb = 0, track, track0 = 0, xtrack, xtrack0 = 0;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <pert.nc> <lon0> <lat0> <overpass.tab>");

  /* Get arguments... */
  lon0 = atof(argv[3]);
  lat0 = atof(argv[4]);

  /* Get control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  orblat = scan_ctl(argc, argv, "ORBLAT", -1, "0", NULL);
  rmax = scan_ctl(argc, argv, "RMAX", -1, "100", NULL);
  obsz = scan_ctl(argc, argv, "OBSZ", -1, "", NULL);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], pertname, pert);

  /* Get Cartesian coordinates... */
  geo2cart(0, lon0, lat0, x0);

  /* Create file... */
  printf("Write overpass data file: %s\n", argv[5]);
  if (!(out = fopen(argv[5], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time (seconds since 2000-01-01T00:00Z)\n"
	  "# $2  = time (string)\n"
	  "# $3  = longitude [deg]\n"
	  "# $4  = latitude [deg]\n"
	  "# $5  = along-track index\n"
	  "# $6  = across-track index\n"
	  "# $7  = orbit number\n"
	  "# $8  = ascending (1=yes, 0=no)\n"
	  "# $9  = scan angle [deg]\n" "# $10 = distance [km]\n\n");

  /* Find nearest footprint... */
  for (track = 0; track < pert->ntrack; track++) {

    /* Check for new orbit... */
    if (track > 0)
      if (pert->lat[track - 1][pert->nxtrack / 2] <= orblat
	  && pert->lat[track][pert->nxtrack / 2] >= orblat) {

	/* Write results... */
	if (sqrt(dmin) <= rmax)
	  write_results(out, pert, track0, xtrack0, orb, dmin, obsz);

	/* Set counters... */
	dmin = 1e100;
	orb++;
      }

    /* Check distance of footprints... */
    for (xtrack = 0; xtrack < pert->nxtrack; xtrack++) {
      geo2cart(0, pert->lon[track][xtrack], pert->lat[track][xtrack], x1);
      if (DIST2(x0, x1) < dmin) {
	dmin = DIST2(x0, x1);
	track0 = track;
	xtrack0 = xtrack;
      }
    }
  }

  /* Write results for last orbit... */
  if (sqrt(dmin) <= rmax)
    write_results(out, pert, track0, xtrack0, orb, dmin, obsz);

  /* Close file... */
  fclose(out);

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void write_results(
  FILE *out,
  pert_t *pert,
  int track0,
  int xtrack0,
  int orb,
  double dmin,
  double obsz) {

  double alpha, xf[3], xs[3], xsf[3], remain;

  int asc, i, year, mon, day, hour, min, sec;

  /* Calculate scan angle... */
  geo2cart(0, pert->lon[track0][xtrack0], pert->lat[track0][xtrack0], xf);
  geo2cart(0, pert->lon[track0][pert->nxtrack / 2],
	   pert->lat[track0][pert->nxtrack / 2], xsf);
  geo2cart(obsz, pert->lon[track0][pert->nxtrack / 2],
	   pert->lat[track0][pert->nxtrack / 2], xs);
  for (i = 0; i < 3; i++) {
    xf[i] -= xs[i];
    xsf[i] -= xs[i];
  }
  alpha = 180. / M_PI * acos(DOTP(xf, xsf) / NORM(xf) / NORM(xsf));
  if (xtrack0 < pert->nxtrack / 2)
    alpha = -alpha;

  /* Get ascending/descending flag... */
  asc = (pert->lat[track0 > 0 ? track0 : track0 + 1][pert->nxtrack / 2]
	 > pert->lat[track0 > 0 ? track0 - 1 : track0][pert->nxtrack / 2]);

  /* Write results... */
  jsec2time(pert->time[track0][xtrack0], &year, &mon, &day,
	    &hour, &min, &sec, &remain);
  fprintf(out,
	  "%.2f %d-%02d-%02dT%02d:%02d:%02dZ %g %g %d %d %d %d %g %g\n",
	  pert->time[track0][xtrack0], year, mon, day, hour, min, sec,
	  pert->lon[track0][xtrack0], pert->lat[track0][xtrack0],
	  track0, xtrack0, orb, asc, alpha, sqrt(dmin));
}
