#include "libairs.h"

/* ------------------------------------------------------------
   Definitions...
   ------------------------------------------------------------ */

/* Maximum number of data points. */
#define NMAX 1000000

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  FILE *in;

  static double x[NMAX], y[NMAX], w[2 * NMAX];

  static size_t n, stride;

  /* Write info... */
  if (argc != 3)
    ERRMSG("Give parameters: <xy.tab> <stride>");

  /* Get stride... */
  stride = (size_t) atoi(argv[2]);

  /* Read data... */
  if (!(in = fopen(argv[1], "r")))
    ERRMSG("Cannot open file!");
  while (fscanf(in, "%lg %lg", &x[n], &y[n]) == 2)
    if (gsl_finite(x[n]) && gsl_finite(y[n]))
      if (++n >= NMAX)
	ERRMSG("Too many data points!");
  fclose(in);

  /* Get statistics... */
  printf("                    n= %lu\n\n", n);
  printf("        xy_covariance= %g\n",
	 gsl_stats_covariance(x, stride, y, stride, n));
  printf("       xy_correlation= %g\n",
	 gsl_stats_correlation(x, stride, y, stride, n));
  printf("          xy_spearman= %g\n",
	 gsl_stats_spearman(x, stride, y, stride, n, w));
  printf("      x_lag1_autocorr= %g\n",
	 gsl_stats_lag1_autocorrelation(x, stride, n));
  printf("      y_lag1_autocorr= %g\n\n",
	 gsl_stats_lag1_autocorrelation(y, stride, n));

  /* Sort data... */
  gsl_sort(x, stride, n);
  gsl_sort(y, stride, n);

  /* Get statistics... */
  printf("               x_mean= %g\n", gsl_stats_mean(x, stride, n));
  printf("              x_sigma= %g\n", gsl_stats_sd(x, stride, n));
  printf("           x_skewness= %g\n", gsl_stats_skew(x, stride, n));
  printf("           x_kurtosis= %g\n\n", gsl_stats_kurtosis(x, stride, n));

  printf("            x_mininum= %g\n", gsl_stats_min(x, stride, n));
  printf("       x_10%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(x, stride, n, 0.1));
  printf("       x_25%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(x, stride, n, 0.25));
  printf("       x_50%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(x, stride, n, 0.5));
  printf("       x_75%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(x, stride, n, 0.75));
  printf("       x_90%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(x, stride, n, 0.9));
  printf("            x_maximum= %g\n\n", gsl_stats_max(x, stride, n));

  printf("        x_absdev_mean= %g\n", gsl_stats_absdev(x, stride, n));
  printf("      x_absdev_median= %g\n",
	 gsl_stats_absdev_m(x, stride, n,
			    gsl_stats_quantile_from_sorted_data(x, stride, n,
								0.5)));
  printf("        x_absdev_zero= %g\n",
	 gsl_stats_absdev_m(x, stride, n, 0.0));
  printf("x_interquartile_range= %g\n\n",
	 gsl_stats_quantile_from_sorted_data(x, stride, n, 0.75)
	 - gsl_stats_quantile_from_sorted_data(x, stride, n, 0.25));

  printf("               y_mean= %g\n", gsl_stats_mean(y, stride, n));
  printf("              y_sigma= %g\n", gsl_stats_sd(y, stride, n));
  printf("           y_skewness= %g\n", gsl_stats_skew(y, stride, n));
  printf("           y_kurtosis= %g\n\n", gsl_stats_kurtosis(y, stride, n));

  printf("            y_mininum= %g\n", gsl_stats_min(y, stride, n));
  printf("       y_10%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(y, stride, n, 0.1));
  printf("       y_25%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(y, stride, n, 0.25));
  printf("       y_50%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(y, stride, n, 0.5));
  printf("       y_75%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(y, stride, n, 0.75));
  printf("       y_90%%_quantile= %g\n",
	 gsl_stats_quantile_from_sorted_data(y, stride, n, 0.9));
  printf("            y_maximum= %g\n\n", gsl_stats_max(y, stride, n));

  printf("        y_absdev_mean= %g\n", gsl_stats_absdev(y, stride, n));
  printf("      y_absdev_median= %g\n",
	 gsl_stats_absdev_m(y, stride, n,
			    gsl_stats_quantile_from_sorted_data(y, stride, n,
								0.5)));
  printf("        y_absdev_zero= %g\n",
	 gsl_stats_absdev_m(y, stride, n, 0.0));
  printf("y_interquartile_range= %g\n\n",
	 gsl_stats_quantile_from_sorted_data(y, stride, n, 0.75)
	 - gsl_stats_quantile_from_sorted_data(y, stride, n, 0.25));


  return (EXIT_SUCCESS);
}
