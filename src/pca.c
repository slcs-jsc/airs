#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;

  static gsl_matrix *a, *v;

  static gsl_vector *s, *w;

  static double lat[AIRS_RAD_GEOTRACK * AIRS_RAD_GEOXTRACK],
    lon[AIRS_RAD_GEOTRACK * AIRS_RAD_GEOXTRACK], mean;

  static size_t channel0, channel1, ichan, itrack, ixtrack, i, j, m, n;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <l1b_file1>");

  /* Get arguments... */
  channel0 = (size_t) scan_ctl(argc, argv, "CHANNEL0", -1, "", NULL);
  channel1 = (size_t) scan_ctl(argc, argv, "CHANNEL1", -1, "", NULL);

  /* Read AIRS data... */
  printf("Read AIRS Level-1B data file: %s\n", argv[2]);
  airs_rad_rdr(argv[2], &airs_rad_gran);

  /* Allocate... */
  m = AIRS_RAD_GEOTRACK * AIRS_RAD_GEOXTRACK;
  n = channel1 - channel0 + 1;
  a = gsl_matrix_calloc(m, n);
  v = gsl_matrix_calloc(n, n);
  s = gsl_vector_calloc(n);
  w = gsl_vector_calloc(n);

  /* Build data matrix... */
  for (itrack = 0; itrack < AIRS_RAD_GEOTRACK; itrack++)
    for (ixtrack = 0; ixtrack < AIRS_RAD_GEOXTRACK; ixtrack++) {
      i = itrack * AIRS_RAD_GEOXTRACK + ixtrack;
      lon[i] = airs_rad_gran.Longitude[itrack][ixtrack];
      lat[i] = airs_rad_gran.Latitude[itrack][ixtrack];
      for (ichan = channel0; ichan <= channel1; ichan++)
	if (airs_rad_gran.radiances[itrack][ixtrack][ichan] > 0)
	  gsl_matrix_set(a, i, (ichan - channel0),
			 brightness(airs_rad_gran.radiances[itrack][ixtrack]
				    [ichan] * 0.001,
				    airs_rad_gran.nominal_freq[ichan]));
    }

  /* Remove column mean... */
  for (j = 0; j < n; j++) {
    mean = 0;
    for (i = 0; i < m; i++)
      mean += gsl_matrix_get(a, i, j) / (double) m;
    printf("mean[%lu] = %g K\n", j, mean);
    for (i = 0; i < m; i++)
      gsl_matrix_set(a, i, j, gsl_matrix_get(a, i, j) - mean);
  }

  /* Calculate SVD... */
  gsl_linalg_SV_decomp(a, v, s, w);

  /*
     https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
   */

  /* Write eigenvalues (variances of PCs)... */
  for (i = 0; i < n; i++)
    printf("lambda_i[%lu] = %g\n", i,
	   gsl_pow_2(gsl_vector_get(s, i)) / ((double) n - 1.0));

  /* Calculate principal components (columns of U x S)... */
  for (j = 0; j < n; j++) {
    printf("\n");
    for (i = 0; i < m; i++)
      printf("%lu %lu %g %g %g %g\n", i, j, lon[i], lat[i],
	     airs_rad_gran.nominal_freq[channel0 + j], gsl_matrix_get(a, i,
								      j) *
	     gsl_vector_get(s, j));
  }

  /* Free... */
  gsl_matrix_free(a);
  gsl_matrix_free(v);
  gsl_vector_free(s);
  gsl_vector_free(w);

  return EXIT_SUCCESS;
}
