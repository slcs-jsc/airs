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
  AIRS Code Collection library definitions.
*/

#include "libairs.h"

/*****************************************************************************/

void add_att(
  int ncid,
  int varid,
  const char *unit,
  const char *long_name) {

  /* Set long name... */
  NC(nc_put_att_text(ncid, varid, "long_name", strlen(long_name), long_name));

  /* Set units... */
  NC(nc_put_att_text(ncid, varid, "units", strlen(unit), unit));
}

/*****************************************************************************/

void add_var(
  int ncid,
  const char *varname,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int *varid,
  int ndims) {

  /* Check if variable exists... */
  if (nc_inq_varid(ncid, varname, varid) != NC_NOERR) {

    /* Define variable... */
    NC(nc_def_var(ncid, varname, type, ndims, dimid, varid));

    /* Set long name... */
    NC(nc_put_att_text
       (ncid, *varid, "long_name", strlen(longname), longname));

    /* Set units... */
    NC(nc_put_att_text(ncid, *varid, "units", strlen(unit), unit));
  }
}

/*****************************************************************************/

void background_poly_help(
  double *xx,
  double *yy,
  int n,
  int dim) {

  double chisq, xx2[WX > WY ? WX : WY], yy2[WX > WY ? WX : WY];

  size_t i, i2, n2 = 0;

  /* Check for nan... */
  for (i = 0; i < (size_t) n; i++)
    if (gsl_finite(yy[i])) {
      xx2[n2] = xx[i];
      yy2[n2] = yy[i];
      n2++;
    }
  if ((int) n2 < dim || n2 < 0.9 * n) {
    for (i = 0; i < (size_t) n; i++)
      yy[i] = GSL_NAN;
    return;
  }

  /* Allocate... */
  gsl_multifit_linear_workspace *work
    = gsl_multifit_linear_alloc((size_t) n2, (size_t) dim);
  gsl_matrix *cov = gsl_matrix_alloc((size_t) dim, (size_t) dim);
  gsl_matrix *X = gsl_matrix_alloc((size_t) n2, (size_t) dim);
  gsl_vector *c = gsl_vector_alloc((size_t) dim);
  gsl_vector *x = gsl_vector_alloc((size_t) n2);
  gsl_vector *y = gsl_vector_alloc((size_t) n2);

  /* Compute polynomial fit... */
  for (i = 0; i < (size_t) n2; i++) {
    gsl_vector_set(x, i, xx2[i]);
    gsl_vector_set(y, i, yy2[i]);
    for (i2 = 0; i2 < (size_t) dim; i2++)
      gsl_matrix_set(X, i, i2, pow(gsl_vector_get(x, i), (double) i2));
  }
  gsl_multifit_linear(X, y, c, cov, &chisq, work);
  for (i = 0; i < (size_t) n; i++)
    yy[i] = gsl_poly_eval(c->data, (int) dim, xx[i]);

  /* Free... */
  gsl_multifit_linear_free(work);
  gsl_matrix_free(cov);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(x);
  gsl_vector_free(y);
}

/*****************************************************************************/

void background_poly(
  wave_t *wave,
  int dim_x,
  int dim_y) {

  double x[WX], x2[WY], y[WX], y2[WY];

  /* Check parameters... */
  if (dim_x <= 0 && dim_y <= 0)
    return;

  /* Copy temperatures to background... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {
      wave->bg[ix][iy] = wave->temp[ix][iy];
      wave->pt[ix][iy] = 0;
    }

  /* Compute fit in x-direction... */
  if (dim_x > 0)
    for (int iy = 0; iy < wave->ny; iy++) {
      for (int ix = 0; ix < wave->nx; ix++) {
	x[ix] = (double) ix;
	y[ix] = wave->bg[ix][iy];
      }
      background_poly_help(x, y, wave->nx, dim_x);
      for (int ix = 0; ix < wave->nx; ix++)
	wave->bg[ix][iy] = y[ix];
    }

  /* Compute fit in y-direction... */
  if (dim_y > 0)
    for (int ix = 0; ix < wave->nx; ix++) {
      for (int iy = 0; iy < wave->ny; iy++) {
	x2[iy] = (int) iy;
	y2[iy] = wave->bg[ix][iy];
      }
      background_poly_help(x2, y2, wave->ny, dim_y);
      for (int iy = 0; iy < wave->ny; iy++)
	wave->bg[ix][iy] = y2[iy];
    }

  /* Recompute perturbations... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++)
      wave->pt[ix][iy] = wave->temp[ix][iy] - wave->bg[ix][iy];
}

/*****************************************************************************/

void background_smooth(
  wave_t *wave,
  int npts_x,
  int npts_y) {

  const double dmax = 2500.0;

  static double help[WX][WY];

  /* Check parameters... */
  if (npts_x <= 0 && npts_y <= 0)
    return;

  /* Smooth background... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {

      /* Init... */
      int n = 0;
      help[ix][iy] = 0;

      /* Set maximum range... */
      const int dx = GSL_MIN(GSL_MIN(npts_x, ix), wave->nx - 1 - ix);
      const int dy = GSL_MIN(GSL_MIN(npts_y, iy), wave->ny - 1 - iy);

      /* Average... */
      for (int i = ix - dx; i <= ix + dx; i++)
	for (int j = iy - dy; j <= iy + dy; j++)
	  if (fabs(wave->x[ix] - wave->x[i]) < dmax &&
	      fabs(wave->y[iy] - wave->y[j]) < dmax) {
	    help[ix][iy] += wave->bg[i][j];
	    n++;
	  }

      /* Normalize... */
      if (n > 0)
	help[ix][iy] /= n;
      else
	help[ix][iy] = GSL_NAN;
    }

  /* Recalculate perturbations... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {
      wave->bg[ix][iy] = help[ix][iy];
      wave->pt[ix][iy] = wave->temp[ix][iy] - wave->bg[ix][iy];
    }
}

/*****************************************************************************/

void create_background(
  wave_t *wave) {

  /* Loop over grid points... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {

      /* Set background for 4.3 micron BT measurements... */
      wave->bg[ix][iy] = 235.626 + 5.38165e-6 * gsl_pow_2(wave->x[ix]
							  -
							  0.5 * (wave->x[0] +
								 wave->x
								 [wave->nx -
								  1]))
	- 1.78519e-12 * gsl_pow_4(wave->x[ix] -
				  0.5 * (wave->x[0] + wave->x[wave->nx - 1]));

      /* Set temperature perturbation... */
      wave->pt[ix][iy] = 0;

      /* Set temperature... */
      wave->temp[ix][iy] = wave->bg[ix][iy];
    }
}

/*****************************************************************************/

void create_noise(
  wave_t *wave,
  double nedt) {

  /* Initialize random number generator... */
  gsl_rng_env_setup();
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, (unsigned long int) time(NULL));

  /* Add noise to temperature... */
  if (nedt > 0)
    for (int ix = 0; ix < wave->nx; ix++)
      for (int iy = 0; iy < wave->ny; iy++)
	wave->temp[ix][iy] += gsl_ran_gaussian(r, nedt);

  /* Free... */
  gsl_rng_free(r);
}

/*****************************************************************************/

void create_wave(
  wave_t *wave,
  double amp,
  double lx,
  double ly,
  double phi,
  double fwhm) {

  /* Loop over grid points... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {

      /* Set wave perturbation... */
      wave->pt[ix][iy] = amp * cos((lx != 0 ? 2 * M_PI / lx : 0) * wave->x[ix]
				   + (ly !=
				      0 ? 2 * M_PI / ly : 0) * wave->y[iy]
				   - phi * M_PI / 180.)
	* (fwhm > 0 ? exp(-0.5 * gsl_pow_2((wave->x[ix]) / (lx * fwhm) * 2.35)
			  -
			  0.5 * gsl_pow_2((wave->y[iy]) / (ly * fwhm) *
					  2.35)) : 1.0);

      /* Add perturbation to temperature... */
      wave->temp[ix][iy] += wave->pt[ix][iy];
    }
}

/*****************************************************************************/

void day2doy(
  int year,
  int mon,
  int day,
  int *doy) {

  const int d0[12] =
    { 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 };
  const int d0l[12] =
    { 1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336 };

  /* Get day of year... */
  if (year % 400 == 0 || (year % 100 != 0 && year % 4 == 0))
    *doy = d0l[mon - 1] + day - 1;
  else
    *doy = d0[mon - 1] + day - 1;
}

/*****************************************************************************/

void doy2day(
  int year,
  int doy,
  int *mon,
  int *day) {

  const int d0[12] =
    { 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 };
  const int d0l[12] =
    { 1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336 };
  int i;

  /* Get month and day... */
  if (year % 400 == 0 || (year % 100 != 0 && year % 4 == 0)) {
    for (i = 11; i >= 0; i--)
      if (d0l[i] <= doy)
	break;
    *mon = i + 1;
    *day = doy - d0l[i] + 1;
  } else {
    for (i = 11; i >= 0; i--)
      if (d0[i] <= doy)
	break;
    *mon = i + 1;
    *day = doy - d0[i] + 1;
  }
}

/*****************************************************************************/

void fit_wave(
  wave_t *wave,
  double amp,
  double phi,
  double kx,
  double ky,
  double *chisq) {

  /* Init... */
  *chisq = 0;

  /* Calculate fit... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {
      wave->fit[ix][iy]
	= amp * cos(kx * wave->x[ix] + ky * wave->y[iy]
		    - phi * M_PI / 180.);
      *chisq += POW2(wave->fit[ix][iy] - wave->pt[ix][iy]);
    }

  /* Calculate chisq... */
  *chisq /= (double) (wave->nx * wave->ny);
}

/*****************************************************************************/

void fft_help(
  double *fcReal,
  double *fcImag,
  int n) {

  double data[2 * PMAX];

  /* Check size... */
  if (n > PMAX)
    ERRMSG("Too many data points!");

  /* Allocate... */
  gsl_fft_complex_wavetable *wavetable
    = gsl_fft_complex_wavetable_alloc((size_t) n);
  gsl_fft_complex_workspace *workspace
    = gsl_fft_complex_workspace_alloc((size_t) n);

  /* Set data (real, complex)... */
  for (int i = 0; i < n; i++) {
    data[2 * i] = fcReal[i];
    data[2 * i + 1] = fcImag[i];
  }

  /* Calculate FFT... */
  gsl_fft_complex_forward(data, 1, (size_t) n, wavetable, workspace);

  /* Copy data... */
  for (int i = 0; i < n; i++) {
    fcReal[i] = data[2 * i];
    fcImag[i] = data[2 * i + 1];
  }

  /* Free... */
  gsl_fft_complex_wavetable_free(wavetable);
  gsl_fft_complex_workspace_free(workspace);
}

/*****************************************************************************/

void fft(
  wave_t *wave,
  double *Amax,
  double *phimax,
  double *lhmax,
  double *kxmax,
  double *kymax,
  double *alphamax,
  double *betamax,
  char *filename) {

  static double A[PMAX][PMAX], phi[PMAX][PMAX], kx[PMAX], ky[PMAX],
    cutReal[PMAX], cutImag[PMAX], boxImag[PMAX][PMAX], boxReal[PMAX][PMAX];

  FILE *out;

  /* Find box... */
  int imin = 9999, jmin = 9999;
  int imax = -9999, jmax = -9999;
  for (int i = 0; i < wave->nx; i++)
    for (int j = 0; j < wave->ny; j++)
      if (gsl_finite(wave->var[i][j])) {
	imin = GSL_MIN(imin, i);
	imax = GSL_MAX(imax, i);
	jmin = GSL_MIN(jmin, j);
	jmax = GSL_MAX(jmax, j);
      }
  const int nx = imax - imin + 1;
  const int ny = jmax - jmin + 1;

  /* Copy data... */
  for (int i = imin; i <= imax; i++)
    for (int j = jmin; j <= jmax; j++) {
      if (gsl_finite(wave->pt[i][j]))
	boxReal[i - imin][j - jmin] = wave->pt[i][j];
      else
	boxReal[i - imin][j - jmin] = 0.0;
      boxImag[i - imin][j - jmin] = 0.0;
    }

  /* FFT of the rows... */
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      cutReal[j] = boxReal[i][j];
      cutImag[j] = boxImag[i][j];
    }
    fft_help(cutReal, cutImag, ny);
    for (int j = 0; j < ny; j++) {
      boxReal[i][j] = cutReal[j];
      boxImag[i][j] = cutImag[j];
    }
  }

  /* FFT of the columns... */
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      cutReal[i] = boxReal[i][j];
      cutImag[i] = boxImag[i][j];
    }
    fft_help(cutReal, cutImag, nx);
    for (int i = 0; i < nx; i++) {
      boxReal[i][j] = cutReal[i];
      boxImag[i][j] = cutImag[i];
    }
  }

  /* Get frequencies, amplitude, and phase... */
  for (int i = 0; i < nx; i++)
    kx[i] = 2. * M_PI * ((i < nx / 2) ? (double) i : -(double) (nx - i))
      / (nx * fabs(wave->x[imax] - wave->x[imin]) / (nx - 1.0));
  for (int j = 0; j < ny; j++)
    ky[j] = 2. * M_PI * ((j < ny / 2) ? (double) j : -(double) (ny - j))
      / (ny * fabs(wave->y[jmax] - wave->y[jmin]) / (ny - 1.0));
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++) {
      A[i][j]
	= (i == 0 && j == 0 ? 1.0 : 2.0) / (nx * ny)
	* sqrt(gsl_pow_2(boxReal[i][j]) + gsl_pow_2(boxImag[i][j]));
      phi[i][j]
	= 180. / M_PI * atan2(boxImag[i][j], boxReal[i][j]);
    }

  /* Check frequencies... */
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      if (kx[i] == 0 || ky[j] == 0) {
	A[i][j] = GSL_NAN;
	phi[i][j] = GSL_NAN;
      }

  /* Find maximum... */
  *Amax = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny / 2; j++)
      if (gsl_finite(A[i][j]) && A[i][j] > *Amax) {
	*Amax = A[i][j];
	*phimax = phi[i][j];
	*kxmax = kx[i];
	*kymax = ky[j];
      }

  /* Get horizontal wavelength... */
  *lhmax = 2 * M_PI / sqrt(gsl_pow_2(*kxmax) + gsl_pow_2(*kymax));

  /* Get propagation direction in xy-plane... */
  *alphamax = 90. - 180. / M_PI * atan2(*kxmax, *kymax);

  /* Get propagation direction in lon,lat-plane... */
  *betamax = *alphamax
    +
    180. / M_PI *
    atan2(wave->lat[wave->nx / 2 >
		    0 ? wave->nx / 2 - 1 : wave->nx / 2][wave->ny / 2]
	  - wave->lat[wave->nx / 2 <
		      wave->nx - 1 ? wave->nx / 2 +
		      1 : wave->nx / 2][wave->ny / 2],
	  wave->lon[wave->nx / 2 >
		    0 ? wave->nx / 2 - 1 : wave->nx / 2][wave->ny / 2]
	  - wave->lon[wave->nx / 2 <
		      wave->nx - 1 ? wave->nx / 2 +
		      1 : wave->nx / 2][wave->ny / 2]);

  /* Save FFT data... */
  if (filename != NULL) {

    /* Write info... */
    printf("Write FFT data: %s\n", filename);

    /* Create file... */
    if (!(out = fopen(filename, "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = altitude [km]\n"
	    "# $2 = wavelength in x-direction [km]\n"
	    "# $3 = wavelength in y-direction [km]\n"
	    "# $4 = wavenumber in x-direction [1/km]\n"
	    "# $5 = wavenumber in y-direction [1/km]\n"
	    "# $6 = amplitude [K]\n" "# $7 = phase [rad]\n");

    /* Write data... */
    for (int i = nx - 1; i > 0; i--) {
      fprintf(out, "\n");
      for (int j = ny / 2; j > 0; j--) {
	int i2 = (i == nx / 2 ? 0 : i);
	int j2 = (j == ny / 2 ? 0 : j);
	fprintf(out, "%g %g %g %g %g %g %g\n", wave->z,
		(kx[i2] != 0 ? 2 * M_PI / kx[i2] : 0),
		(ky[j2] != 0 ? 2 * M_PI / ky[j2] : 0),
		kx[i2], ky[j2], A[i2][j2], phi[i2][j2]);
      }
    }

    /* Close file... */
    fclose(out);
  }
}

/*****************************************************************************/

void gauss(
  wave_t *wave,
  double fwhm) {

  static double help[WX][WY];

  /* Check parameters... */
  if (fwhm <= 0)
    return;

  /* Compute sigma^2... */
  const double sigma2 = gsl_pow_2(fwhm / 2.3548);

  /* Loop over data points... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {

      /* Init... */
      double wsum = 0;
      help[ix][iy] = 0;

      /* Average... */
      for (int ix2 = 0; ix2 < wave->nx; ix2++)
	for (int iy2 = 0; iy2 < wave->ny; iy2++) {
	  const double d2 = gsl_pow_2(wave->x[ix] - wave->x[ix2])
	    + gsl_pow_2(wave->y[iy] - wave->y[iy2]);
	  if (d2 <= 9 * sigma2) {
	    const double w = exp(-d2 / (2 * sigma2));
	    wsum += w;
	    help[ix][iy] += w * wave->pt[ix2][iy2];
	  }
	}

      /* Normalize... */
      wave->pt[ix][iy] = help[ix][iy] / wsum;
    }
}

/*****************************************************************************/

void hamming(
  wave_t *wave,
  int niter) {

  static double help[WX][WY];

  /* Iterations... */
  for (int iter = 0; iter < niter; iter++) {

    /* Filter in x direction... */
    for (int ix = 0; ix < wave->nx; ix++)
      for (int iy = 0; iy < wave->ny; iy++)
	help[ix][iy]
	  = 0.23 * wave->pt[ix > 0 ? ix - 1 : ix][iy]
	  + 0.54 * wave->pt[ix][iy]
	  + 0.23 * wave->pt[ix < wave->nx - 1 ? ix + 1 : ix][iy];

    /* Filter in y direction... */
    for (int ix = 0; ix < wave->nx; ix++)
      for (int iy = 0; iy < wave->ny; iy++)
	wave->pt[ix][iy]
	  = 0.23 * help[ix][iy > 0 ? iy - 1 : iy]
	  + 0.54 * help[ix][iy]
	  + 0.23 * help[ix][iy < wave->ny - 1 ? iy + 1 : iy];
  }
}

/*****************************************************************************/

void intpol_x(
  wave_t *wave,
  int n) {

  double dummy, x[WX], xc[WX][3], xc2[WX][3], y[WX];

  /* Check parameters... */
  if (n <= 0)
    return;
  if (n > WX)
    ERRMSG("Too many data points!");

  /* Set new x-coordinates... */
  for (int i = 0; i < n; i++)
    x[i] = LIN(0.0, wave->x[0], n - 1.0, wave->x[wave->nx - 1], i);

  /* Allocate... */
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline
    = gsl_spline_alloc(gsl_interp_cspline, (size_t) wave->nx);

  /* Loop over scans... */
  for (int iy = 0; iy < wave->ny; iy++) {

    /* Interpolate Cartesian coordinates... */
    for (int ix = 0; ix < wave->nx; ix++)
      geo2cart(0, wave->lon[ix][iy], wave->lat[ix][iy], xc[ix]);
    for (int ic = 0; ic < 3; ic++) {
      for (int ix = 0; ix < wave->nx; ix++)
	y[ix] = xc[ix][ic];
      gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
      for (int i = 0; i < n; i++)
	xc2[i][ic] = gsl_spline_eval(spline, x[i], acc);
    }
    for (int i = 0; i < n; i++)
      cart2geo(xc2[i], &dummy, &wave->lon[i][iy], &wave->lat[i][iy]);

    /* Interpolate temperature... */
    for (int ix = 0; ix < wave->nx; ix++)
      y[ix] = wave->temp[ix][iy];
    gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
    for (int i = 0; i < n; i++)
      wave->temp[i][iy] = gsl_spline_eval(spline, x[i], acc);

    /* Interpolate background... */
    for (int ix = 0; ix < wave->nx; ix++)
      y[ix] = wave->bg[ix][iy];
    gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
    for (int i = 0; i < n; i++)
      wave->bg[i][iy] = gsl_spline_eval(spline, x[i], acc);

    /* Interpolate perturbations... */
    for (int ix = 0; ix < wave->nx; ix++)
      y[ix] = wave->pt[ix][iy];
    gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
    for (int i = 0; i < n; i++)
      wave->pt[i][iy] = gsl_spline_eval(spline, x[i], acc);

    /* Interpolate variance... */
    for (int ix = 0; ix < wave->nx; ix++)
      y[ix] = wave->var[ix][iy];
    gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
    for (int i = 0; i < n; i++)
      wave->var[i][iy] = gsl_spline_eval(spline, x[i], acc);
  }

  /* Free... */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  /* Set new x-coordinates... */
  for (int i = 0; i < n; i++)
    wave->x[i] = x[i];
  wave->nx = n;
}

/*****************************************************************************/

void median(
  wave_t *wave,
  int dx) {

  static double data[WX * WY], help[WX][WY];

  /* Check parameters... */
  if (dx <= 0)
    return;

  /* Loop over data points... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {

      /* Init... */
      size_t n = 0;

      /* Get data... */
      for (int ix2 = GSL_MAX(ix - dx, 0);
	   ix2 < GSL_MIN(ix + dx, wave->nx - 1); ix2++)
	for (int iy2 = GSL_MAX(iy - dx, 0);
	     iy2 < GSL_MIN(iy + dx, wave->ny - 1); iy2++) {
	  data[n] = wave->pt[ix2][iy2];
	  n++;
	}

      /* Normalize... */
      gsl_sort(data, 1, n);
      help[ix][iy] = gsl_stats_median_from_sorted_data(data, 1, n);
    }

  /* Loop over data points... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++)
      wave->pt[ix][iy] = help[ix][iy];
}

/*****************************************************************************/

void merge_y(
  wave_t *wave1,
  wave_t *wave2) {

  /* Check data... */
  if (wave1->nx != wave2->nx)
    ERRMSG("Across-track sizes do not match!");
  if (wave1->ny + wave2->ny > WY)
    ERRMSG("Too many data points!");

  /* Get offset in y direction... */
  double y =
    wave1->y[wave1->ny - 1] + (wave1->y[wave1->ny - 1] -
			       wave1->y[0]) / (wave1->ny - 1);

  /* Merge data... */
  for (int ix = 0; ix < wave2->nx; ix++)
    for (int iy = 0; iy < wave2->ny; iy++) {
      wave1->y[wave1->ny + iy] = y + wave2->y[iy];
      wave1->lon[ix][wave1->ny + iy] = wave2->lon[ix][iy];
      wave1->lat[ix][wave1->ny + iy] = wave2->lat[ix][iy];
      wave1->temp[ix][wave1->ny + iy] = wave2->temp[ix][iy];
      wave1->bg[ix][wave1->ny + iy] = wave2->bg[ix][iy];
      wave1->pt[ix][wave1->ny + iy] = wave2->pt[ix][iy];
      wave1->var[ix][wave1->ny + iy] = wave2->var[ix][iy];
    }

  /* Increment counter... */
  wave1->ny += wave2->ny;
}

/*****************************************************************************/

void noise(
  wave_t *wave,
  double *mu,
  double *sig) {

  /* Init... */
  int n = 0;
  *mu = 0;
  *sig = 0;

  /* Estimate noise (Immerkaer, 1996)... */
  for (int ix = 1; ix < wave->nx - 1; ix++)
    for (int iy = 1; iy < wave->ny - 1; iy++) {

      /* Check data... */
      int okay = 1;
      for (int ix2 = ix - 1; ix2 <= ix + 1; ix2++)
	for (int iy2 = iy - 1; iy2 <= iy + 1; iy2++)
	  if (!gsl_finite(wave->temp[ix2][iy2]))
	    okay = 0;
      if (!okay)
	continue;

      /* Get mean noise... */
      n++;
      *mu += wave->temp[ix][iy];
      *sig += gsl_pow_2(+4. / 6. * wave->temp[ix][iy]
			- 2. / 6. * (wave->temp[ix - 1][iy]
				     + wave->temp[ix + 1][iy]
				     + wave->temp[ix][iy - 1]
				     + wave->temp[ix][iy + 1])
			+ 1. / 6. * (wave->temp[ix - 1][iy - 1]
				     + wave->temp[ix + 1][iy - 1]
				     + wave->temp[ix - 1][iy + 1]
				     + wave->temp[ix + 1][iy + 1]));
    }

  /* Normalize... */
  *mu /= (double) n;
  *sig = sqrt(*sig / (double) n);
}

/*****************************************************************************/

void period(
  wave_t *wave,
  double lxymax,
  double dlxy,
  double *Amax,
  double *phimax,
  double *lhmax,
  double *kxmax,
  double *kymax,
  double *alphamax,
  double *betamax,
  char *filename) {

  FILE *out;

  static double kx[PMAX], ky[PMAX], A[PMAX][PMAX], phi[PMAX][PMAX],
    cx[PMAX][WX], cy[PMAX][WY], sx[PMAX][WX], sy[PMAX][WY];

  int imin, imax, jmin, jmax, lmax = 0, mmax = 0;

  /* Compute wavenumbers and periodogram coefficients... */
  for (double lx = -lxymax; lx <= lxymax; lx += dlxy) {
    kx[lmax] = (lx != 0 ? 2 * M_PI / lx : 0);
    for (int i = 0; i < wave->nx; i++) {
      cx[lmax][i] = cos(kx[lmax] * wave->x[i]);
      sx[lmax][i] = sin(kx[lmax] * wave->x[i]);
    }
    if ((++lmax) > PMAX)
      ERRMSG("Too many wavenumbers for periodogram!");
  }
  for (double ly = 0; ly <= lxymax; ly += dlxy) {
    ky[mmax] = (ly != 0 ? 2 * M_PI / ly : 0);
    for (int j = 0; j < wave->ny; j++) {
      cy[mmax][j] = cos(ky[mmax] * wave->y[j]);
      sy[mmax][j] = sin(ky[mmax] * wave->y[j]);
    }
    if ((++mmax) > PMAX)
      ERRMSG("Too many wavenumbers for periodogram!");
  }

  /* Find area... */
  imin = jmin = 9999;
  imax = jmax = -9999;
  for (int i = 0; i < wave->nx; i++)
    for (int j = 0; j < wave->ny; j++)
      if (gsl_finite(wave->var[i][j])) {
	imin = GSL_MIN(imin, i);
	imax = GSL_MAX(imax, i);
	jmin = GSL_MIN(jmin, j);
	jmax = GSL_MAX(jmax, j);
      }

  /* Get Nyquist frequencies... */
  const double kx_ny =
    M_PI / fabs((wave->x[imax] - wave->x[imin]) /
		((double) imax - (double) imin));
  const double ky_ny =
    M_PI / fabs((wave->y[jmax] - wave->y[jmin]) /
		((double) jmax - (double) jmin));

  /* Loop over wavelengths... */
  for (int l = 0; l < lmax; l++)
    for (int m = 0; m < mmax; m++) {

      /* Check frequencies... */
      if (kx[l] == 0 || fabs(kx[l]) > kx_ny ||
	  ky[m] == 0 || fabs(ky[m]) > ky_ny) {
	A[l][m] = GSL_NAN;
	phi[l][m] = GSL_NAN;
	continue;
      }

      /* Compute periodogram... */
      double a = 0, b = 0, c = 0;
      for (int i = imin; i <= imax; i++)
	for (int j = jmin; j <= jmax; j++)
	  if (gsl_finite(wave->var[i][j])) {
	    a += wave->pt[i][j] * (cx[l][i] * cy[m][j] - sx[l][i] * sy[m][j]);
	    b += wave->pt[i][j] * (sx[l][i] * cy[m][j] + cx[l][i] * sy[m][j]);
	    c++;
	  }
      a *= 2. / c;
      b *= 2. / c;

      /* Get amplitude and phase... */
      A[l][m] = sqrt(gsl_pow_2(a) + gsl_pow_2(b));
      phi[l][m] = atan2(b, a) * 180. / M_PI;
    }

  /* Find maximum... */
  *Amax = 0;
  for (int l = 0; l < lmax; l++)
    for (int m = 0; m < mmax; m++)
      if (gsl_finite(A[l][m]) && A[l][m] > *Amax) {
	*Amax = A[l][m];
	*phimax = phi[l][m];
	*kxmax = kx[l];
	*kymax = ky[m];
      }

  /* Get horizontal wavelength... */
  *lhmax = 2 * M_PI / sqrt(gsl_pow_2(*kxmax) + gsl_pow_2(*kymax));

  /* Get propagation direction in xy-plane... */
  *alphamax = 90. - 180. / M_PI * atan2(*kxmax, *kymax);

  /* Get propagation direction in lon,lat-plane... */
  *betamax = *alphamax
    +
    180. / M_PI *
    atan2(wave->lat[wave->nx / 2 >
		    0 ? wave->nx / 2 - 1 : wave->nx / 2][wave->ny / 2]
	  - wave->lat[wave->nx / 2 <
		      wave->nx - 1 ? wave->nx / 2 +
		      1 : wave->nx / 2][wave->ny / 2],
	  wave->lon[wave->nx / 2 >
		    0 ? wave->nx / 2 - 1 : wave->nx / 2][wave->ny / 2]
	  - wave->lon[wave->nx / 2 <
		      wave->nx - 1 ? wave->nx / 2 +
		      1 : wave->nx / 2][wave->ny / 2]);

  /* Save periodogram data... */
  if (filename != NULL) {

    /* Write info... */
    printf("Write periodogram data: %s\n", filename);

    /* Create file... */
    if (!(out = fopen(filename, "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = altitude [km]\n"
	    "# $2 = wavelength in x-direction [km]\n"
	    "# $3 = wavelength in y-direction [km]\n"
	    "# $4 = wavenumber in x-direction [1/km]\n"
	    "# $5 = wavenumber in y-direction [1/km]\n"
	    "# $6 = amplitude [K]\n" "# $7 = phase [rad]\n");

    /* Write data... */
    for (int l = 0; l < lmax; l++) {
      fprintf(out, "\n");
      for (int m = 0; m < mmax; m++)
	fprintf(out, "%g %g %g %g %g %g %g\n", wave->z,
		(kx[l] != 0 ? 2 * M_PI / kx[l] : 0),
		(ky[m] != 0 ? 2 * M_PI / ky[m] : 0),
		kx[l], ky[m], A[l][m], phi[l][m]);
    }

    /* Close file... */
    fclose(out);
  }
}

/*****************************************************************************/

void pert2wave(
  pert_t *pert,
  wave_t *wave,
  int track0,
  int track1,
  int xtrack0,
  int xtrack1) {

  double x0[3], x1[3];

  /* Check ranges... */
  track0 = GSL_MIN(GSL_MAX(track0, 0), pert->ntrack - 1);
  track1 = GSL_MIN(GSL_MAX(track1, 0), pert->ntrack - 1);
  xtrack0 = GSL_MIN(GSL_MAX(xtrack0, 0), pert->nxtrack - 1);
  xtrack1 = GSL_MIN(GSL_MAX(xtrack1, 0), pert->nxtrack - 1);

  /* Set size... */
  wave->nx = xtrack1 - xtrack0 + 1;
  if (wave->nx > WX)
    ERRMSG("Too many across-track values!");
  wave->ny = track1 - track0 + 1;
  if (wave->ny > WY)
    ERRMSG("Too many along-track values!");

  /* Loop over footprints... */
  for (int itrack = track0; itrack <= track1; itrack++)
    for (int ixtrack = xtrack0; ixtrack <= xtrack1; ixtrack++) {

      /* Get distances... */
      if (itrack == track0) {
	wave->x[0] = 0;
	if (ixtrack > xtrack0) {
	  geo2cart(0, pert->lon[itrack][ixtrack - 1],
		   pert->lat[itrack][ixtrack - 1], x0);
	  geo2cart(0, pert->lon[itrack][ixtrack],
		   pert->lat[itrack][ixtrack], x1);
	  wave->x[ixtrack - xtrack0] =
	    wave->x[ixtrack - xtrack0 - 1] + DIST(x0, x1);
	}
      }
      if (ixtrack == xtrack0) {
	wave->y[0] = 0;
	if (itrack > track0) {
	  geo2cart(0, pert->lon[itrack - 1][ixtrack],
		   pert->lat[itrack - 1][ixtrack], x0);
	  geo2cart(0, pert->lon[itrack][ixtrack],
		   pert->lat[itrack][ixtrack], x1);
	  wave->y[itrack - track0] =
	    wave->y[itrack - track0 - 1] + DIST(x0, x1);
	}
      }

      /* Save geolocation... */
      wave->time = pert->time[(track0 + track1) / 2][(xtrack0 + xtrack1) / 2];
      wave->z = 0;
      wave->lon[ixtrack - xtrack0][itrack - track0] =
	pert->lon[itrack][ixtrack];
      wave->lat[ixtrack - xtrack0][itrack - track0] =
	pert->lat[itrack][ixtrack];

      /* Save temperature data... */
      wave->temp[ixtrack - xtrack0][itrack - track0]
	= pert->bt[itrack][ixtrack];
      wave->bg[ixtrack - xtrack0][itrack - track0]
	= pert->bt[itrack][ixtrack] - pert->pt[itrack][ixtrack];
      wave->pt[ixtrack - xtrack0][itrack - track0]
	= pert->pt[itrack][ixtrack];
      wave->var[ixtrack - xtrack0][itrack - track0]
	= pert->var[itrack][ixtrack];
    }
}

/*****************************************************************************/

void read_l1(
  char *filename,
  airs_l1_t *l1) {

  int ncid, varid;

  /* Open netCDF file... */
  printf("Read AIRS Level-1 file: %s\n", filename);
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Read data... */
  NC(nc_inq_varid(ncid, "l1_time", &varid));
  NC(nc_get_var_double(ncid, varid, l1->time[0]));
  NC(nc_inq_varid(ncid, "l1_lon", &varid));
  NC(nc_get_var_double(ncid, varid, l1->lon[0]));
  NC(nc_inq_varid(ncid, "l1_lat", &varid));
  NC(nc_get_var_double(ncid, varid, l1->lat[0]));
  NC(nc_inq_varid(ncid, "l1_sat_z", &varid));
  NC(nc_get_var_double(ncid, varid, l1->sat_z));
  NC(nc_inq_varid(ncid, "l1_sat_lon", &varid));
  NC(nc_get_var_double(ncid, varid, l1->sat_lon));
  NC(nc_inq_varid(ncid, "l1_sat_lat", &varid));
  NC(nc_get_var_double(ncid, varid, l1->sat_lat));
  NC(nc_inq_varid(ncid, "l1_nu", &varid));
  NC(nc_get_var_double(ncid, varid, l1->nu));
  NC(nc_inq_varid(ncid, "l1_rad", &varid));
  NC(nc_get_var_float(ncid, varid, l1->rad[0][0]));

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_l2(
  char *filename,
  airs_l2_t *l2) {

  int ncid, varid;

  /* Open netCDF file... */
  printf("Read AIRS Level-2 file: %s\n", filename);
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Read data... */
  NC(nc_inq_varid(ncid, "l2_time", &varid));
  NC(nc_get_var_double(ncid, varid, l2->time[0]));
  NC(nc_inq_varid(ncid, "l2_z", &varid));
  NC(nc_get_var_double(ncid, varid, l2->z[0][0]));
  NC(nc_inq_varid(ncid, "l2_lon", &varid));
  NC(nc_get_var_double(ncid, varid, l2->lon[0]));
  NC(nc_inq_varid(ncid, "l2_lat", &varid));
  NC(nc_get_var_double(ncid, varid, l2->lat[0]));
  NC(nc_inq_varid(ncid, "l2_press", &varid));
  NC(nc_get_var_double(ncid, varid, l2->p));
  NC(nc_inq_varid(ncid, "l2_temp", &varid));
  NC(nc_get_var_double(ncid, varid, l2->t[0][0]));

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_pert(
  char *filename,
  char *pertname,
  pert_t *pert) {

  static char varname[LEN];

  static int dimid[2], ncid, varid;

  static size_t ntrack, nxtrack, start[2] = { 0, 0 }, count[2] = { 1, 1 };

  /* Write info... */
  printf("Read perturbation data: %s\n", filename);

  /* Open netCDF file... */
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "NTRACK", &dimid[0]));
  NC(nc_inq_dimid(ncid, "NXTRACK", &dimid[1]));
  NC(nc_inq_dimlen(ncid, dimid[0], &ntrack));
  NC(nc_inq_dimlen(ncid, dimid[1], &nxtrack));
  if (nxtrack > PERT_NXTRACK)
    ERRMSG("Too many tracks!");
  if (ntrack > PERT_NTRACK)
    ERRMSG("Too many scans!");
  pert->ntrack = (int) ntrack;
  pert->nxtrack = (int) nxtrack;
  count[1] = nxtrack;

  /* Read data... */
  NC(nc_inq_varid(ncid, "time", &varid));
  for (size_t itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->time[itrack]));
  }

  NC(nc_inq_varid(ncid, "lon", &varid));
  for (size_t itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->lon[itrack]));
  }

  NC(nc_inq_varid(ncid, "lat", &varid));
  for (size_t itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->lat[itrack]));
  }

  NC(nc_inq_varid(ncid, "bt_8mu", &varid));
  for (size_t itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->dc[itrack]));
  }

  sprintf(varname, "bt_%s", pertname);
  NC(nc_inq_varid(ncid, varname, &varid));
  for (size_t itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->bt[itrack]));
  }

  sprintf(varname, "bt_%s_pt", pertname);
  NC(nc_inq_varid(ncid, varname, &varid));
  for (size_t itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->pt[itrack]));
  }

  sprintf(varname, "bt_%s_var", pertname);
  NC(nc_inq_varid(ncid, varname, &varid));
  for (size_t itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->var[itrack]));
  }

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_retr(
  char *filename,
  ret_t *ret) {

  static double help[NDS * NPG];

  int dimid, ids = 0, ncid, varid;

  size_t nds, np, ntrack, nxtrack;

  /* Write info... */
  printf("Read retrieval data: %s\n", filename);

  /* Open netCDF file... */
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Read new retrieval file format... */
  if (nc_inq_dimid(ncid, "L1_NTRACK", &dimid) == NC_NOERR) {

    /* Get dimensions... */
    NC(nc_inq_dimid(ncid, "RET_NP", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &np));
    ret->np = (int) np;
    if (ret->np > NPG)
      ERRMSG("Too many data points!");

    NC(nc_inq_dimid(ncid, "L1_NTRACK", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &ntrack));
    NC(nc_inq_dimid(ncid, "L1_NXTRACK", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &nxtrack));
    ret->nds = (int) (ntrack * nxtrack);
    if (ret->nds > NDS)
      ERRMSG("Too many data sets!");

    /* Read time... */
    NC(nc_inq_varid(ncid, "l1_time", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    ids = 0;
    for (size_t itrack = 0; itrack < ntrack; itrack++)
      for (size_t ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (int ip = 0; ip < ret->np; ip++)
	  ret->time[ids][ip] = help[ids];
	ids++;
      }

    /* Read altitudes... */
    NC(nc_inq_varid(ncid, "ret_z", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    ids = 0;
    for (size_t itrack = 0; itrack < ntrack; itrack++)
      for (size_t ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (int ip = 0; ip < ret->np; ip++)
	  ret->z[ids][ip] = help[ip];
	ids++;
      }

    /* Read longitudes... */
    NC(nc_inq_varid(ncid, "l1_lon", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    ids = 0;
    for (size_t itrack = 0; itrack < ntrack; itrack++)
      for (size_t ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (int ip = 0; ip < ret->np; ip++)
	  ret->lon[ids][ip] = help[ids];
	ids++;
      }

    /* Read latitudes... */
    NC(nc_inq_varid(ncid, "l1_lat", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    ids = 0;
    for (size_t itrack = 0; itrack < ntrack; itrack++)
      for (size_t ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (int ip = 0; ip < ret->np; ip++)
	  ret->lat[ids][ip] = help[ids];
	ids++;
      }

    /* Read temperatures... */
    NC(nc_inq_varid(ncid, "ret_temp", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    ids = 0;
    for (size_t itrack = 0; itrack < ntrack; itrack++)
      for (size_t ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (int ip = 0; ip < ret->np; ip++)
	  ret->t[ids][ip] =
	    help[(itrack * nxtrack + ixtrack) * (size_t) np + (size_t) ip];
	ids++;
      }
  }

  /* Read old retrieval file format... */
  if (nc_inq_dimid(ncid, "np", &dimid) == NC_NOERR) {

    /* Get dimensions... */
    NC(nc_inq_dimid(ncid, "np", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &np));
    ret->np = (int) np;
    if (ret->np > NPG)
      ERRMSG("Too many data points!");

    NC(nc_inq_dimid(ncid, "nds", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &nds));
    ret->nds = (int) nds;
    if (ret->nds > NDS)
      ERRMSG("Too many data sets!");

    /* Read data... */
    NC(nc_inq_varid(ncid, "time", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->time);

    NC(nc_inq_varid(ncid, "z", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->z);

    NC(nc_inq_varid(ncid, "lon", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->lon);

    NC(nc_inq_varid(ncid, "lat", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->lat);

    NC(nc_inq_varid(ncid, "press", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->p);

    NC(nc_inq_varid(ncid, "temp", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->t);

    NC(nc_inq_varid(ncid, "temp_apr", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->t_apr);

    NC(nc_inq_varid(ncid, "temp_total", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->t_tot);

    NC(nc_inq_varid(ncid, "temp_noise", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->t_noise);

    NC(nc_inq_varid(ncid, "temp_formod", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->t_fm);

    NC(nc_inq_varid(ncid, "temp_cont", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->t_cont);

    NC(nc_inq_varid(ncid, "temp_res", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    read_retr_help(help, ret->nds, ret->np, ret->t_res);

    NC(nc_inq_varid(ncid, "chisq", &varid));
    NC(nc_get_var_double(ncid, varid, ret->chisq));
  }

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_retr_help(
  double *help,
  int nds,
  int np,
  double mat[NDS][NPG]) {

  int n = 0;

  for (int ip = 0; ip < np; ip++)
    for (int ids = 0; ids < nds; ids++)
      mat[ids][ip] = help[n++];
}

/*****************************************************************************/

void read_wave(
  char *filename,
  wave_t *wave) {

  FILE *in;

  char line[LEN];

  double rtime, rz, rlon, rlat, rx, ry, ryold = -1e10,
    rtemp, rbg, rpt, rvar, rfit;

  /* Init... */
  wave->nx = 0;
  wave->ny = 0;

  /* Write info... */
  printf("Read wave data: %s\n", filename);

  /* Open file... */
  if (!(in = fopen(filename, "r")))
    ERRMSG("Cannot open file!");

  /* Read data... */
  while (fgets(line, LEN, in))
    if (sscanf(line, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &rtime,
	       &rz, &rlon, &rlat, &rx, &ry, &rtemp, &rbg, &rpt,
	       &rvar, &rfit) == 11) {

      /* Set index... */
      if (ry != ryold) {
	if ((++wave->ny >= WY))
	  ERRMSG("Too many y-values!");
	wave->nx = 0;
      } else if ((++wave->nx) >= WX)
	ERRMSG("Too many x-values!");
      ryold = ry;

      /* Save data... */
      wave->time = rtime;
      wave->z = rz;
      wave->lon[wave->nx][wave->ny] = rlon;
      wave->lat[wave->nx][wave->ny] = rlat;
      wave->x[wave->nx] = rx;
      wave->y[wave->ny] = ry;
      wave->temp[wave->nx][wave->ny] = rtemp;
      wave->bg[wave->nx][wave->ny] = rbg;
      wave->pt[wave->nx][wave->ny] = rpt;
      wave->var[wave->nx][wave->ny] = rvar;
      wave->var[wave->nx][wave->ny] = rfit;
    }

  /* Increment counters... */
  wave->nx++;
  wave->ny++;

  /* Close file... */
  fclose(in);
}

/*****************************************************************************/

void rad2wave(
  airs_rad_gran_t *gran,
  double *nu,
  int nd,
  wave_t *wave) {

  double x0[3], x1[3];

  int ichan[AIRS_RAD_CHANNEL];

  /* Get channel numbers... */
  for (int id = 0; id < nd; id++) {
    for (ichan[id] = 0; ichan[id] < AIRS_RAD_CHANNEL; ichan[id]++)
      if (fabs(gran->nominal_freq[ichan[id]] - nu[id]) < 0.1)
	break;
    if (ichan[id] >= AIRS_RAD_CHANNEL)
      ERRMSG("Could not find channel!");
  }

  /* Set size... */
  wave->nx = AIRS_RAD_GEOXTRACK;
  wave->ny = AIRS_RAD_GEOTRACK;
  if (wave->nx > WX || wave->ny > WY)
    ERRMSG("Wave struct too small!");

  /* Set Cartesian coordinates... */
  geo2cart(0, gran->Longitude[0][0], gran->Latitude[0][0], x0);
  for (int xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
    geo2cart(0, gran->Longitude[0][xtrack], gran->Latitude[0][xtrack], x1);
    wave->x[xtrack] = DIST(x0, x1);
  }
  for (int track = 0; track < AIRS_RAD_GEOTRACK; track++) {
    geo2cart(0, gran->Longitude[track][0], gran->Latitude[track][0], x1);
    wave->y[track] = DIST(x0, x1);
  }

  /* Set geolocation... */
  wave->time =
    gran->Time[AIRS_RAD_GEOTRACK / 2][AIRS_RAD_GEOXTRACK / 2] - 220838400;
  wave->z = 0;
  for (int track = 0; track < AIRS_RAD_GEOTRACK; track++)
    for (int xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
      wave->lon[xtrack][track] = gran->Longitude[track][xtrack];
      wave->lat[xtrack][track] = gran->Latitude[track][xtrack];
    }

  /* Set brightness temperature... */
  for (int track = 0; track < AIRS_RAD_GEOTRACK; track++)
    for (int xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
      wave->temp[xtrack][track] = 0;
      wave->bg[xtrack][track] = 0;
      wave->pt[xtrack][track] = 0;
      wave->var[xtrack][track] = 0;
      for (int id = 0; id < nd; id++) {
	if ((gran->state[track][xtrack] != 0)
	    || (gran->ExcludedChans[ichan[id]] > 2)
	    || (gran->CalChanSummary[ichan[id]] & 8)
	    || (gran->CalChanSummary[ichan[id]] & (32 + 64))
	    || (gran->CalFlag[track][ichan[id]] & 16))
	  wave->temp[xtrack][track] = GSL_NAN;
	else
	  wave->temp[xtrack][track]
	    += BRIGHT(gran->radiances[track][xtrack][ichan[id]] * 1e-3,
		      gran->nominal_freq[ichan[id]]) / nd;
      }
    }
}

/*****************************************************************************/

void ret2wave(
  ret_t *ret,
  wave_t *wave,
  int dataset,
  int ip) {

  /* Initialize... */
  wave->nx = 90;
  if (wave->nx > WX)
    ERRMSG("Too many across-track values!");
  wave->ny = 135;
  if (wave->ny > WY)
    ERRMSG("Too many along-track values!");
  if (ip < 0 || ip >= ret->np)
    ERRMSG("Altitude index out of range!");

  /* Loop over data sets and data points... */
  for (int ids = 0; ids < ret->nds; ids++) {

    /* Get horizontal indices... */
    const int ix = ids % 90;
    const int iy = ids / 90;

    /* Get distances... */
    double x0[3], x1[3];
    if (iy == 0) {
      geo2cart(0.0, ret->lon[0][0], ret->lat[0][0], x0);
      geo2cart(0.0, ret->lon[ids][ip], ret->lat[ids][ip], x1);
      wave->x[ix] = DIST(x0, x1);
    }
    if (ix == 0) {
      geo2cart(0.0, ret->lon[0][0], ret->lat[0][0], x0);
      geo2cart(0.0, ret->lon[ids][ip], ret->lat[ids][ip], x1);
      wave->y[iy] = DIST(x0, x1);
    }

    /* Save geolocation... */
    wave->time = ret->time[0][0];
    if (ix == 0 && iy == 0)
      wave->z = ret->z[ids][ip];
    wave->lon[ix][iy] = ret->lon[ids][ip];
    wave->lat[ix][iy] = ret->lat[ids][ip];

    /* Save temperature... */
    if (dataset == 1)
      wave->temp[ix][iy] = ret->t[ids][ip];
    else if (dataset == 2)
      wave->temp[ix][iy] = ret->t_apr[ids][ip];
  }
}

/*****************************************************************************/

void variance(
  wave_t *wave,
  double dh) {

  /* Check parameters... */
  if (dh <= 0)
    return;

  /* Compute squared radius... */
  const double dh2 = gsl_pow_2(dh);

  /* Get sampling distances... */
  const int dx =
    (int) (dh / fabs(wave->x[wave->nx - 1] - wave->x[0]) * (wave->nx - 1.0) +
	   1);
  const int dy =
    (int) (dh / fabs(wave->y[wave->ny - 1] - wave->y[0]) * (wave->ny - 1.0) +
	   1);

  /* Loop over data points... */
  for (int ix = 0; ix < wave->nx; ix++)
    for (int iy = 0; iy < wave->ny; iy++) {

      /* Init... */
      double mu = 0, help = 0;
      int n = 0;

      /* Get data... */
      for (int ix2 = GSL_MAX(ix - dx, 0);
	   ix2 <= GSL_MIN(ix + dx, wave->nx - 1); ix2++)
	for (int iy2 = GSL_MAX(iy - dy, 0);
	     iy2 <= GSL_MIN(iy + dy, wave->ny - 1); iy2++)
	  if ((gsl_pow_2(wave->x[ix] - wave->x[ix2])
	       + gsl_pow_2(wave->y[iy] - wave->y[iy2])) <= dh2)
	    if (gsl_finite(wave->pt[ix2][iy2])) {
	      mu += wave->pt[ix2][iy2];
	      help += gsl_pow_2(wave->pt[ix2][iy2]);
	      n++;
	    }

      /* Compute local variance... */
      if (n > 1)
	wave->var[ix][iy] = help / n - gsl_pow_2(mu / n);
      else
	wave->var[ix][iy] = GSL_NAN;
    }
}

/*****************************************************************************/

void write_l1(
  char *filename,
  airs_l1_t *l1) {

  int dimid[10], ncid, time_id, lon_id, lat_id,
    sat_z_id, sat_lon_id, sat_lat_id, nu_id, rad_id;

  /* Open or create netCDF file... */
  printf("Write AIRS Level-1 file: %s\n", filename);
  if (nc_open(filename, NC_WRITE, &ncid) != NC_NOERR) {
    NC(nc_create(filename, NC_CLOBBER, &ncid));
  } else {
    NC(nc_redef(ncid));
  }

  /* Set dimensions... */
  if (nc_inq_dimid(ncid, "L1_NTRACK", &dimid[0]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L1_NTRACK", L1_NTRACK, &dimid[0]));
  if (nc_inq_dimid(ncid, "L1_NXTRACK", &dimid[1]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L1_NXTRACK", L1_NXTRACK, &dimid[1]));
  if (nc_inq_dimid(ncid, "L1_NCHAN", &dimid[2]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L1_NCHAN", L1_NCHAN, &dimid[2]));

  /* Add variables... */
  add_var(ncid, "l1_time", "s", "time (seconds since 2000-01-01T00:00Z)",
	  NC_DOUBLE, dimid, &time_id, 2);
  add_var(ncid, "l1_lon", "deg", "longitude", NC_DOUBLE, dimid, &lon_id, 2);
  add_var(ncid, "l1_lat", "deg", "latitude", NC_DOUBLE, dimid, &lat_id, 2);
  add_var(ncid, "l1_sat_z", "km", "satellite altitude",
	  NC_DOUBLE, dimid, &sat_z_id, 1);
  add_var(ncid, "l1_sat_lon", "deg", "satellite longitude",
	  NC_DOUBLE, dimid, &sat_lon_id, 1);
  add_var(ncid, "l1_sat_lat", "deg", "satellite latitude",
	  NC_DOUBLE, dimid, &sat_lat_id, 1);
  add_var(ncid, "l1_nu", "cm^-1", "channel wavenumber",
	  NC_DOUBLE, &dimid[2], &nu_id, 1);
  add_var(ncid, "l1_rad", "W/(m^2 sr cm^-1)", "channel radiance",
	  NC_FLOAT, dimid, &rad_id, 3);

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Write data... */
  NC(nc_put_var_double(ncid, time_id, l1->time[0]));
  NC(nc_put_var_double(ncid, lon_id, l1->lon[0]));
  NC(nc_put_var_double(ncid, lat_id, l1->lat[0]));
  NC(nc_put_var_double(ncid, sat_z_id, l1->sat_z));
  NC(nc_put_var_double(ncid, sat_lon_id, l1->sat_lon));
  NC(nc_put_var_double(ncid, sat_lat_id, l1->sat_lat));
  NC(nc_put_var_double(ncid, nu_id, l1->nu));
  NC(nc_put_var_float(ncid, rad_id, l1->rad[0][0]));

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void write_l2(
  char *filename,
  airs_l2_t *l2) {

  int dimid[10], ncid, time_id, z_id, lon_id, lat_id, p_id, t_id;

  /* Create netCDF file... */
  printf("Write AIRS Level-2 file: %s\n", filename);
  if (nc_open(filename, NC_WRITE, &ncid) != NC_NOERR) {
    NC(nc_create(filename, NC_CLOBBER, &ncid));
  } else {
    NC(nc_redef(ncid));
  }

  /* Set dimensions... */
  if (nc_inq_dimid(ncid, "L2_NTRACK", &dimid[0]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L2_NTRACK", L2_NTRACK, &dimid[0]));
  if (nc_inq_dimid(ncid, "L2_NXTRACK", &dimid[1]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L2_NXTRACK", L2_NXTRACK, &dimid[1]));
  if (nc_inq_dimid(ncid, "L2_NLAY", &dimid[2]) != NC_NOERR)
    NC(nc_def_dim(ncid, "L2_NLAY", L2_NLAY, &dimid[2]));

  /* Add variables... */
  add_var(ncid, "l2_time", "s", "time (seconds since 2000-01-01T00:00Z)",
	  NC_DOUBLE, dimid, &time_id, 2);
  add_var(ncid, "l2_z", "km", "altitude", NC_DOUBLE, dimid, &z_id, 3);
  add_var(ncid, "l2_lon", "deg", "longitude", NC_DOUBLE, dimid, &lon_id, 2);
  add_var(ncid, "l2_lat", "deg", "latitude", NC_DOUBLE, dimid, &lat_id, 2);
  add_var(ncid, "l2_press", "hPa", "pressure",
	  NC_DOUBLE, &dimid[2], &p_id, 1);
  add_var(ncid, "l2_temp", "K", "temperature", NC_DOUBLE, dimid, &t_id, 3);

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Write data... */
  NC(nc_put_var_double(ncid, time_id, l2->time[0]));
  NC(nc_put_var_double(ncid, z_id, l2->z[0][0]));
  NC(nc_put_var_double(ncid, lon_id, l2->lon[0]));
  NC(nc_put_var_double(ncid, lat_id, l2->lat[0]));
  NC(nc_put_var_double(ncid, p_id, l2->p));
  NC(nc_put_var_double(ncid, t_id, l2->t[0][0]));

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void write_wave(
  char *filename,
  wave_t *wave) {

  FILE *out;

  /* Write info... */
  printf("Write wave data: %s\n", filename);

  /* Create file... */
  if (!(out = fopen(filename, "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time (seconds since 2000-01-01T00:00Z)\n"
	  "# $2  = altitude [km]\n"
	  "# $3  = longitude [deg]\n"
	  "# $4  = latitude [deg]\n"
	  "# $5  = across-track distance [km]\n"
	  "# $6  = along-track distance [km]\n"
	  "# $7  = temperature [K]\n"
	  "# $8  = background [K]\n"
	  "# $9  = perturbation [K]\n"
	  "# $10 = variance [K^2]\n" "# $11 = fitting model [K]\n");

  /* Write data... */
  for (int j = 0; j < wave->ny; j++) {
    fprintf(out, "\n");
    for (int i = 0; i < wave->nx; i++)
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g\n",
	      wave->time, wave->z, wave->lon[i][j], wave->lat[i][j],
	      wave->x[i], wave->y[j], wave->temp[i][j], wave->bg[i][j],
	      wave->pt[i][j], wave->var[i][j], wave->fit[i][j]);
  }

  /* Close file... */
  fclose(out);
}
