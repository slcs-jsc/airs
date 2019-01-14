#include "libairs.h"

/* ------------------------------------------------------------
   Definitions...
   ------------------------------------------------------------ */

/* Maximum number of data points. */
#define NMAX 1000000

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  FILE *in;

  static double x[NMAX], y[NMAX], y2[NMAX], y3[NMAX],
    t, work[2 * NMAX], v0, t0, t1;

  static int i, idx, it, n, n2, nt, rm;

  /* Write info... */
  if (argc != 6)
    ERRMSG("Give parameters: <xy.tab> <rm> var0> <t0> <t1>");
  rm = atoi(argv[2]);
  v0 = atof(argv[3]);
  t0 = atof(argv[4]);
  t1 = atof(argv[5]);

  /* Read data... */
  printf("# Read time series: %s\n\n", argv[1]);
  if (!(in = fopen(argv[1], "r")))
    ERRMSG("Cannot open file!");
  while (fscanf(in, "%lg %lg", &x[n], &y[n]) == 2)
    if (x[n] >= t0 && x[n] <= t1)
      if (++n >= NMAX)
	ERRMSG("Too many data points!");
  fclose(in);

  /* Interpolate time series... */
  for (t = x[0]; t <= x[n - 1]; t += 86400) {
    idx = locate(x, n, t);
    y2[n2] = LIN(x[idx], y[idx], x[idx + 1], y[idx + 1], t);
    if (++n2 >= NMAX)
      ERRMSG("Too many data points!");
  }

  /* Remove running mean... */
  for (i = 0; i < n2; i++) {
    y3[i] = 0;
    nt = 0;
    for (it = -rm; it <= rm; it++)
      if (i + it >= 0 && i + it < n2) {
	y3[i] += y2[i + it];
	nt++;
      }
    y3[i] /= nt;
  }
  nt = n2;
  n2 = 0;
  for (i = 0; i < nt; i++)
    if (y3[i] > v0) {
      y2[n2] = y2[i] - y3[i];
      n2++;
    }

  /* Loop over time lag... */
  for (it = 0; it < n2; it++) {

    /* Shift time series... */
    for (i = 0; i < n2; i++)
      y3[i] = y2[i + it < n2 ? i + it : i + it - n2];

    /* Get correlation coefficient... */
    printf("%d %g %g\n", it,
	   gsl_stats_correlation(y2, 1, y3, 1, (size_t) n2),
	   gsl_stats_spearman(y2, 1, y3, 1, (size_t) n2, work));
  }

  return (EXIT_SUCCESS);
}
