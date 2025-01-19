/*
  This file is part of the AIRS Code Collection.
  
  the AIRS Code Collections is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.
  
  The AIRS Code Collection is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with the AIRS Code Collection. If not, see
  <http://www.gnu.org/licenses/>.
  
  Copyright (C) 2019-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Estimate horizontal wavelength sensitivity.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static double chisq;

  /* Check arguments... */
  if (argc != 8)
    ERRMSG("Give parameters: <width> <n> <lxmin> <lxmax> <dlx> <fwhm> <dim>");

  /* Get arguments... */
  const double width = atof(argv[1]);
  const int n = atoi(argv[2]);
  const double lxmin = atof(argv[3]);
  const double lxmax = atof(argv[4]);
  const double dlx = atoi(argv[5]);
  const double fwhm = atof(argv[6]);
  const int dim = atoi(argv[7]);

  /* Initialize... */
  gsl_vector *c = gsl_vector_alloc((size_t) dim);
  gsl_matrix *cov = gsl_matrix_alloc((size_t) dim, (size_t) dim);
  gsl_multifit_linear_workspace *work
    = gsl_multifit_linear_alloc((size_t) n, (size_t) dim);
  gsl_matrix *X = gsl_matrix_alloc((size_t) n, (size_t) dim);
  gsl_vector *xvec = gsl_vector_alloc((size_t) n);
  gsl_vector *yvec = gsl_vector_alloc((size_t) n);
  gsl_vector *yfit = gsl_vector_alloc((size_t) n);

  /* Loop over wavelengths... */
  for (double lx = lxmin; lx <= lxmax; lx += dlx) {

    /* Initialize... */
    double vmean = 0;
    double vmean2 = 0;

    /* Loop over phases... */
    for (double phi = 0; phi < 2 * M_PI; phi += M_PI / 180) {

      /* Initialize... */
      double var = 0, var2 = 0, wsum = 0;

      /* Set wave... */
      for (int i = 0; i < n; i++) {
	gsl_vector_set(xvec, (size_t) i, width / (n - 1.0) * i - width / 2.);
	gsl_vector_set(yvec, (size_t) i,
		       sin(2 * M_PI / lx * gsl_vector_get(xvec, (size_t) i) +
			   phi));
	if (fwhm > 0) {
	  double w = gsl_ran_gaussian_pdf(gsl_vector_get(xvec, (size_t) i),
					  fwhm * lx / 2.3548);
	  gsl_vector_set(yvec, (size_t) i,
			 w * gsl_vector_get(yvec, (size_t) i));
	  wsum += w;
	}
      }
      if (wsum > 0)
	gsl_vector_scale(yvec, 1 / wsum);

      /* Detrending... */
      for (int i = 0; i < n; i++)
	for (int i2 = 0; i2 < dim; i2++)
	  gsl_matrix_set(X, (size_t) i, (size_t) i2,
			 pow(gsl_vector_get(xvec, (size_t) i), 1. * i2));
      gsl_multifit_linear(X, yvec, c, cov, &chisq, work);
      for (int i = 0; i < n; i++)
	gsl_vector_set(yfit, (size_t) i, gsl_vector_get(yvec, (size_t) i)
		       - gsl_poly_eval(c->data, (int) dim,
				       gsl_vector_get(xvec, (size_t) i)));

      /* Compute variances... */
      for (int i = 0; i < n; i++) {
	var += gsl_pow_2(gsl_vector_get(yfit, (size_t) i)) / (double) n;
	var2 += gsl_pow_2(gsl_vector_get(yvec, (size_t) i)) / (double) n;
      }
      vmean += var;
      vmean2 += var2;
    }

    /* Write output... */
    printf("%g %g\n", lx, 100 * vmean / vmean2);
  }

  return EXIT_SUCCESS;
}
