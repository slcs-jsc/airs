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

  gsl_multifit_linear_workspace *work;
  gsl_matrix *cov, *X;
  gsl_vector *c, *x, *y;

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
  work = gsl_multifit_linear_alloc((size_t) n2, (size_t) dim);
  cov = gsl_matrix_alloc((size_t) dim, (size_t) dim);
  X = gsl_matrix_alloc((size_t) n2, (size_t) dim);
  c = gsl_vector_alloc((size_t) dim);
  x = gsl_vector_alloc((size_t) n2);
  y = gsl_vector_alloc((size_t) n2);

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
  wave_t * wave,
  int dim_x,
  int dim_y) {

  double x[WX], x2[WY], y[WX], y2[WY];

  int ix, iy;

  /* Copy temperatures to background... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {
      wave->bg[ix][iy] = wave->temp[ix][iy];
      wave->pt[ix][iy] = 0;
    }

  /* Check parameters... */
  if (dim_x <= 0 && dim_y <= 0)
    return;

  /* Compute fit in x-direction... */
  if (dim_x > 0)
    for (iy = 0; iy < wave->ny; iy++) {
      for (ix = 0; ix < wave->nx; ix++) {
	x[ix] = (double) ix;
	y[ix] = wave->bg[ix][iy];
      }
      background_poly_help(x, y, wave->nx, dim_x);
      for (ix = 0; ix < wave->nx; ix++)
	wave->bg[ix][iy] = y[ix];
    }

  /* Compute fit in y-direction... */
  if (dim_y > 0)
    for (ix = 0; ix < wave->nx; ix++) {
      for (iy = 0; iy < wave->ny; iy++) {
	x2[iy] = (int) iy;
	y2[iy] = wave->bg[ix][iy];
      }
      background_poly_help(x2, y2, wave->ny, dim_y);
      for (iy = 0; iy < wave->ny; iy++)
	wave->bg[ix][iy] = y2[iy];
    }

  /* Recompute perturbations... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++)
      wave->pt[ix][iy] = wave->temp[ix][iy] - wave->bg[ix][iy];
}

/*****************************************************************************/

void background_smooth(
  wave_t * wave,
  int npts_x,
  int npts_y) {

  static double help[WX][WY], dmax = 2500.;

  int dx, dy, i, j, ix, iy, n;

  /* Check parameters... */
  if (npts_x <= 0 && npts_y <= 0)
    return;

  /* Smooth background... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {

      /* Init... */
      n = 0;
      help[ix][iy] = 0;

      /* Set maximum range... */
      dx = GSL_MIN(GSL_MIN(npts_x, ix), wave->nx - 1 - ix);
      dy = GSL_MIN(GSL_MIN(npts_y, iy), wave->ny - 1 - iy);

      /* Average... */
      for (i = ix - dx; i <= ix + dx; i++)
	for (j = iy - dy; j <= iy + dy; j++)
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
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {
      wave->bg[ix][iy] = help[ix][iy];
      wave->pt[ix][iy] = wave->temp[ix][iy] - wave->bg[ix][iy];
    }
}

/*****************************************************************************/

void create_background(
  wave_t * wave) {

  int ix, iy;

  /* Loop over grid points... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {

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
  wave_t * wave,
  double nedt) {

  gsl_rng *r;

  int ix, iy;

  /* Initialize random number generator... */
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, (unsigned long int) time(NULL));

  /* Add noise to temperature... */
  if (nedt > 0)
    for (ix = 0; ix < wave->nx; ix++)
      for (iy = 0; iy < wave->ny; iy++)
	wave->temp[ix][iy] += gsl_ran_gaussian(r, nedt);

  /* Free... */
  gsl_rng_free(r);
}

/*****************************************************************************/

void create_wave(
  wave_t * wave,
  double amp,
  double lx,
  double ly,
  double phi,
  double fwhm) {

  int ix, iy;

  /* Loop over grid points... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {

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

  int d0[12] = { 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 };
  int d0l[12] = { 1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336 };

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

  int d0[12] = { 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 };
  int d0l[12] = { 1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336 };
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

void fft_help(
  double *fcReal,
  double *fcImag,
  int n) {

  gsl_fft_complex_wavetable *wavetable;
  gsl_fft_complex_workspace *workspace;

  double data[2 * PMAX];

  int i;

  /* Check size... */
  if (n > PMAX)
    ERRMSG("Too many data points!");

  /* Allocate... */
  wavetable = gsl_fft_complex_wavetable_alloc((size_t) n);
  workspace = gsl_fft_complex_workspace_alloc((size_t) n);

  /* Set data (real, complex)... */
  for (i = 0; i < n; i++) {
    data[2 * i] = fcReal[i];
    data[2 * i + 1] = fcImag[i];
  }

  /* Calculate FFT... */
  gsl_fft_complex_forward(data, 1, (size_t) n, wavetable, workspace);

  /* Copy data... */
  for (i = 0; i < n; i++) {
    fcReal[i] = data[2 * i];
    fcImag[i] = data[2 * i + 1];
  }

  /* Free... */
  gsl_fft_complex_wavetable_free(wavetable);
  gsl_fft_complex_workspace_free(workspace);
}

/*****************************************************************************/

void fft(
  wave_t * wave,
  double *Amax,
  double *phimax,
  double *lhmax,
  double *alphamax,
  double *betamax,
  char *filename) {

  static double A[PMAX][PMAX], phi[PMAX][PMAX], kx[PMAX], ky[PMAX],
    kxmax, kymax, cutReal[PMAX], cutImag[PMAX],
    boxImag[PMAX][PMAX], boxReal[PMAX][PMAX];

  FILE *out;

  int i, i2, imin, imax, j, j2, jmin, jmax, nx, ny;

  /* Find box... */
  imin = jmin = 9999;
  imax = jmax = -9999;
  for (i = 0; i < wave->nx; i++)
    for (j = 0; j < wave->ny; j++)
      if (gsl_finite(wave->var[i][j])) {
	imin = GSL_MIN(imin, i);
	imax = GSL_MAX(imax, i);
	jmin = GSL_MIN(jmin, j);
	jmax = GSL_MAX(jmax, j);
      }
  nx = imax - imin + 1;
  ny = jmax - jmin + 1;

  /* Copy data... */
  for (i = imin; i <= imax; i++)
    for (j = jmin; j <= jmax; j++) {
      if (gsl_finite(wave->pt[i][j]))
	boxReal[i - imin][j - jmin] = wave->pt[i][j];
      else
	boxReal[i - imin][j - jmin] = 0.0;
      boxImag[i - imin][j - jmin] = 0.0;
    }

  /* FFT of the rows... */
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      cutReal[j] = boxReal[i][j];
      cutImag[j] = boxImag[i][j];
    }
    fft_help(cutReal, cutImag, ny);
    for (j = 0; j < ny; j++) {
      boxReal[i][j] = cutReal[j];
      boxImag[i][j] = cutImag[j];
    }
  }

  /* FFT of the columns... */
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      cutReal[i] = boxReal[i][j];
      cutImag[i] = boxImag[i][j];
    }
    fft_help(cutReal, cutImag, nx);
    for (i = 0; i < nx; i++) {
      boxReal[i][j] = cutReal[i];
      boxImag[i][j] = cutImag[i];
    }
  }

  /* Get frequencies, amplitude, and phase... */
  for (i = 0; i < nx; i++)
    kx[i] = 2. * M_PI * ((i < nx / 2) ? (double) i : -(double) (nx - i))
      / (nx * fabs(wave->x[imax] - wave->x[imin]) / (nx - 1.0));
  for (j = 0; j < ny; j++)
    ky[j] = 2. * M_PI * ((j < ny / 2) ? (double) j : -(double) (ny - j))
      / (ny * fabs(wave->y[jmax] - wave->y[jmin]) / (ny - 1.0));
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      A[i][j]
	= (i == 0 && j == 0 ? 1.0 : 2.0) / (nx * ny)
	* sqrt(gsl_pow_2(boxReal[i][j]) + gsl_pow_2(boxImag[i][j]));
      phi[i][j]
	= 180. / M_PI * atan2(boxImag[i][j], boxReal[i][j]);
    }

  /* Check frequencies... */
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      if (kx[i] == 0 || ky[j] == 0) {
	A[i][j] = GSL_NAN;
	phi[i][j] = GSL_NAN;
      }

  /* Find maximum... */
  *Amax = 0;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny / 2; j++)
      if (gsl_finite(A[i][j]) && A[i][j] > *Amax) {
	*Amax = A[i][j];
	*phimax = phi[i][j];
	kxmax = kx[i];
	kymax = ky[j];
	imax = i;
	jmax = j;
      }

  /* Get horizontal wavelength... */
  *lhmax = 2 * M_PI / sqrt(gsl_pow_2(kxmax) + gsl_pow_2(kymax));

  /* Get propagation direction in xy-plane... */
  *alphamax = 90. - 180. / M_PI * atan2(kxmax, kymax);

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
    for (i = nx - 1; i > 0; i--) {
      fprintf(out, "\n");
      for (j = ny / 2; j > 0; j--) {
	i2 = (i == nx / 2 ? 0 : i);
	j2 = (j == ny / 2 ? 0 : j);
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
  wave_t * wave,
  double fwhm) {

  static double d2, help[WX][WY], sigma2, w, wsum;

  int ix, ix2, iy, iy2;

  /* Check parameters... */
  if (fwhm <= 0)
    return;

  /* Compute sigma^2... */
  sigma2 = gsl_pow_2(fwhm / 2.3548);

  /* Loop over data points... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {

      /* Init... */
      wsum = 0;
      help[ix][iy] = 0;

      /* Average... */
      for (ix2 = 0; ix2 < wave->nx; ix2++)
	for (iy2 = 0; iy2 < wave->ny; iy2++) {
	  d2 = gsl_pow_2(wave->x[ix] - wave->x[ix2])
	    + gsl_pow_2(wave->y[iy] - wave->y[iy2]);
	  if (d2 <= 9 * sigma2) {
	    w = exp(-d2 / (2 * sigma2));
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
  wave_t * wave,
  int niter) {

  static double help[WX][WY];

  int iter, ix, iy;

  /* Iterations... */
  for (iter = 0; iter < niter; iter++) {

    /* Filter in x direction... */
    for (ix = 0; ix < wave->nx; ix++)
      for (iy = 0; iy < wave->ny; iy++)
	help[ix][iy]
	  = 0.23 * wave->pt[ix > 0 ? ix - 1 : ix][iy]
	  + 0.54 * wave->pt[ix][iy]
	  + 0.23 * wave->pt[ix < wave->nx - 1 ? ix + 1 : ix][iy];

    /* Filter in y direction... */
    for (ix = 0; ix < wave->nx; ix++)
      for (iy = 0; iy < wave->ny; iy++)
	wave->pt[ix][iy]
	  = 0.23 * help[ix][iy > 0 ? iy - 1 : iy]
	  + 0.54 * help[ix][iy]
	  + 0.23 * help[ix][iy < wave->ny - 1 ? iy + 1 : iy];
  }
}

/*****************************************************************************/

void intpol_x(
  wave_t * wave,
  int n) {

  gsl_interp_accel *acc;
  gsl_spline *spline;

  double dummy, x[WX], xc[WX][3], xc2[WX][3], y[WX];

  int i, ic, ix, iy;

  /* Check parameters... */
  if (n <= 0)
    return;
  if (n > WX)
    ERRMSG("Too many data points!");

  /* Set new x-coordinates... */
  for (i = 0; i < n; i++)
    x[i] = LIN(0.0, wave->x[0], n - 1.0, wave->x[wave->nx - 1], i);

  /* Allocate... */
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, (size_t) wave->nx);

  /* Loop over scans... */
  for (iy = 0; iy < wave->ny; iy++) {

    /* Interpolate Cartesian coordinates... */
    for (ix = 0; ix < wave->nx; ix++)
      geo2cart(0, wave->lon[ix][iy], wave->lat[ix][iy], xc[ix]);
    for (ic = 0; ic < 3; ic++) {
      for (ix = 0; ix < wave->nx; ix++)
	y[ix] = xc[ix][ic];
      gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
      for (i = 0; i < n; i++)
	xc2[i][ic] = gsl_spline_eval(spline, x[i], acc);
    }
    for (i = 0; i < n; i++)
      cart2geo(xc2[i], &dummy, &wave->lon[i][iy], &wave->lat[i][iy]);

    /* Interpolate temperature... */
    for (ix = 0; ix < wave->nx; ix++)
      y[ix] = wave->temp[ix][iy];
    gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
    for (i = 0; i < n; i++)
      wave->temp[i][iy] = gsl_spline_eval(spline, x[i], acc);

    /* Interpolate background... */
    for (ix = 0; ix < wave->nx; ix++)
      y[ix] = wave->bg[ix][iy];
    gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
    for (i = 0; i < n; i++)
      wave->bg[i][iy] = gsl_spline_eval(spline, x[i], acc);

    /* Interpolate perturbations... */
    for (ix = 0; ix < wave->nx; ix++)
      y[ix] = wave->pt[ix][iy];
    gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
    for (i = 0; i < n; i++)
      wave->pt[i][iy] = gsl_spline_eval(spline, x[i], acc);

    /* Interpolate variance... */
    for (ix = 0; ix < wave->nx; ix++)
      y[ix] = wave->var[ix][iy];
    gsl_spline_init(spline, wave->x, y, (size_t) wave->nx);
    for (i = 0; i < n; i++)
      wave->var[i][iy] = gsl_spline_eval(spline, x[i], acc);
  }

  /* Free... */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  /* Set new x-coordinates... */
  for (i = 0; i < n; i++)
    wave->x[i] = x[i];
  wave->nx = n;
}

/*****************************************************************************/

void median(
  wave_t * wave,
  int dx) {

  static double data[WX * WY], help[WX][WY];

  int ix, ix2, iy, iy2;

  size_t n;

  /* Check parameters... */
  if (dx <= 0)
    return;

  /* Loop over data points... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {

      /* Init... */
      n = 0;

      /* Get data... */
      for (ix2 = GSL_MAX(ix - dx, 0); ix2 < GSL_MIN(ix + dx, wave->nx - 1);
	   ix2++)
	for (iy2 = GSL_MAX(iy - dx, 0); iy2 < GSL_MIN(iy + dx, wave->ny - 1);
	     iy2++) {
	  data[n] = wave->pt[ix2][iy2];
	  n++;
	}

      /* Normalize... */
      gsl_sort(data, 1, n);
      help[ix][iy] = gsl_stats_median_from_sorted_data(data, 1, n);
    }

  /* Loop over data points... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++)
      wave->pt[ix][iy] = help[ix][iy];
}

/*****************************************************************************/

void merge_y(
  wave_t * wave1,
  wave_t * wave2) {

  double y;

  int ix, iy;

  /* Check data... */
  if (wave1->nx != wave2->nx)
    ERRMSG("Across-track sizes do not match!");
  if (wave1->ny + wave2->ny > WY)
    ERRMSG("Too many data points!");

  /* Get offset in y direction... */
  y =
    wave1->y[wave1->ny - 1] + (wave1->y[wave1->ny - 1] -
			       wave1->y[0]) / (wave1->ny - 1);

  /* Merge data... */
  for (ix = 0; ix < wave2->nx; ix++)
    for (iy = 0; iy < wave2->ny; iy++) {
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
  wave_t * wave,
  double *mu,
  double *sig) {

  int ix, ix2, iy, iy2, n = 0, okay;

  /* Init... */
  *mu = 0;
  *sig = 0;

  /* Estimate noise (Immerkaer, 1996)... */
  for (ix = 1; ix < wave->nx - 1; ix++)
    for (iy = 1; iy < wave->ny - 1; iy++) {

      /* Check data... */
      okay = 1;
      for (ix2 = ix - 1; ix2 <= ix + 1; ix2++)
	for (iy2 = iy - 1; iy2 <= iy + 1; iy2++)
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
  wave_t * wave,
  double *Amax,
  double *phimax,
  double *lhmax,
  double *alphamax,
  double *betamax,
  char *filename) {

  FILE *out;

  static double kx[PMAX], ky[PMAX], kx_ny, ky_ny, kxmax, kymax, A[PMAX][PMAX],
    phi[PMAX][PMAX], cx[PMAX][WX], cy[PMAX][WY], sx[PMAX][WX], sy[PMAX][WY],
    a, b, c, lx, ly, lxymax = 1000, dlxy = 10;

  int i, imin, imax, j, jmin, jmax, l, lmax = 0, m, mmax = 0;

  /* Compute wavenumbers and periodogram coefficients... */
  for (lx = -lxymax; lx <= lxymax; lx += dlxy) {
    kx[lmax] = (lx != 0 ? 2 * M_PI / lx : 0);
    for (i = 0; i < wave->nx; i++) {
      cx[lmax][i] = cos(kx[lmax] * wave->x[i]);
      sx[lmax][i] = sin(kx[lmax] * wave->x[i]);
    }
    if ((++lmax) > PMAX)
      ERRMSG("Too many wavenumbers for periodogram!");
  }
  for (ly = 0; ly <= lxymax; ly += dlxy) {
    ky[mmax] = (ly != 0 ? 2 * M_PI / ly : 0);
    for (j = 0; j < wave->ny; j++) {
      cy[mmax][j] = cos(ky[mmax] * wave->y[j]);
      sy[mmax][j] = sin(ky[mmax] * wave->y[j]);
    }
    if ((++mmax) > PMAX)
      ERRMSG("Too many wavenumbers for periodogram!");
  }

  /* Find area... */
  imin = jmin = 9999;
  imax = jmax = -9999;
  for (i = 0; i < wave->nx; i++)
    for (j = 0; j < wave->ny; j++)
      if (gsl_finite(wave->var[i][j])) {
	imin = GSL_MIN(imin, i);
	imax = GSL_MAX(imax, i);
	jmin = GSL_MIN(jmin, j);
	jmax = GSL_MAX(jmax, j);
      }

  /* Get Nyquist frequencies... */
  kx_ny =
    M_PI / fabs((wave->x[imax] - wave->x[imin]) /
		((double) imax - (double) imin));
  ky_ny =
    M_PI / fabs((wave->y[jmax] - wave->y[jmin]) /
		((double) jmax - (double) jmin));

  /* Loop over wavelengths... */
  for (l = 0; l < lmax; l++)
    for (m = 0; m < mmax; m++) {

      /* Check frequencies... */
      if (kx[l] == 0 || fabs(kx[l]) > kx_ny ||
	  ky[m] == 0 || fabs(ky[m]) > ky_ny) {
	A[l][m] = GSL_NAN;
	phi[l][m] = GSL_NAN;
	continue;
      }

      /* Compute periodogram... */
      a = b = c = 0;
      for (i = imin; i <= imax; i++)
	for (j = jmin; j <= jmax; j++)
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
  for (l = 0; l < lmax; l++)
    for (m = 0; m < mmax; m++)
      if (gsl_finite(A[l][m]) && A[l][m] > *Amax) {
	*Amax = A[l][m];
	*phimax = phi[l][m];
	kxmax = kx[l];
	kymax = ky[m];
	imax = i;
	jmax = j;
      }

  /* Get horizontal wavelength... */
  *lhmax = 2 * M_PI / sqrt(gsl_pow_2(kxmax) + gsl_pow_2(kymax));

  /* Get propagation direction in xy-plane... */
  *alphamax = 90. - 180. / M_PI * atan2(kxmax, kymax);

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
    for (l = 0; l < lmax; l++) {
      fprintf(out, "\n");
      for (m = 0; m < mmax; m++)
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
  pert_t * pert,
  wave_t * wave,
  int track0,
  int track1,
  int xtrack0,
  int xtrack1) {

  double x0[3], x1[3];

  int itrack, ixtrack;

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
  for (itrack = track0; itrack <= track1; itrack++)
    for (ixtrack = xtrack0; ixtrack <= xtrack1; ixtrack++) {

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
  airs_l1_t * l1) {

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
  airs_l2_t * l2) {

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
  pert_t * pert) {

  static char varname[LEN];

  static int dimid[2], ncid, varid;

  static size_t itrack, ntrack, nxtrack, start[2] = { 0, 0 }, count[2] = {
    1, 1
  };

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
  for (itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->time[itrack]));
  }

  NC(nc_inq_varid(ncid, "lon", &varid));
  for (itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->lon[itrack]));
  }

  NC(nc_inq_varid(ncid, "lat", &varid));
  for (itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->lat[itrack]));
  }

  NC(nc_inq_varid(ncid, "bt_8mu", &varid));
  for (itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->dc[itrack]));
  }

  sprintf(varname, "bt_%s", pertname);
  NC(nc_inq_varid(ncid, varname, &varid));
  for (itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->bt[itrack]));
  }

  sprintf(varname, "bt_%s_pt", pertname);
  NC(nc_inq_varid(ncid, varname, &varid));
  for (itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->pt[itrack]));
  }

  sprintf(varname, "bt_%s_var", pertname);
  NC(nc_inq_varid(ncid, varname, &varid));
  for (itrack = 0; itrack < ntrack; itrack++) {
    start[0] = itrack;
    NC(nc_get_vara_double(ncid, varid, start, count, pert->var[itrack]));
  }

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_retr(
  char *filename,
  ret_t * ret) {

  static double help[NDS * NPG];

  int dimid, ids = 0, ip, ncid, varid;

  size_t itrack, ixtrack, nds, np, ntrack, nxtrack;

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
    for (itrack = 0; itrack < ntrack; itrack++)
      for (ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (ip = 0; ip < ret->np; ip++)
	  ret->time[ids][ip] = help[ids];
	ids++;
      }

    /* Read altitudes... */
    NC(nc_inq_varid(ncid, "ret_z", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    ids = 0;
    for (itrack = 0; itrack < ntrack; itrack++)
      for (ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (ip = 0; ip < ret->np; ip++)
	  ret->z[ids][ip] = help[ip];
	ids++;
      }

    /* Read longitudes... */
    NC(nc_inq_varid(ncid, "l1_lon", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    ids = 0;
    for (itrack = 0; itrack < ntrack; itrack++)
      for (ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (ip = 0; ip < ret->np; ip++)
	  ret->lon[ids][ip] = help[ids];
	ids++;
      }

    /* Read latitudes... */
    NC(nc_inq_varid(ncid, "l1_lat", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    ids = 0;
    for (itrack = 0; itrack < ntrack; itrack++)
      for (ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (ip = 0; ip < ret->np; ip++)
	  ret->lat[ids][ip] = help[ids];
	ids++;
      }

    /* Read temperatures... */
    NC(nc_inq_varid(ncid, "ret_temp", &varid));
    NC(nc_get_var_double(ncid, varid, help));
    ids = 0;
    for (itrack = 0; itrack < ntrack; itrack++)
      for (ixtrack = 0; ixtrack < nxtrack; ixtrack++) {
	for (ip = 0; ip < ret->np; ip++)
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

  int ids, ip, n = 0;

  for (ip = 0; ip < np; ip++)
    for (ids = 0; ids < nds; ids++)
      mat[ids][ip] = help[n++];
}

/*****************************************************************************/

void read_wave(
  char *filename,
  wave_t * wave) {

  FILE *in;

  char line[LEN];

  double rtime, rz, rlon, rlat, rx, ry, ryold = -1e10, rtemp, rbg, rpt, rvar;

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
    if (sscanf(line, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &rtime,
	       &rz, &rlon, &rlat, &rx, &ry, &rtemp, &rbg, &rpt,
	       &rvar) == 10) {

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
    }

  /* Increment counters... */
  wave->nx++;
  wave->ny++;

  /* Close file... */
  fclose(in);
}

/*****************************************************************************/

void rad2wave(
  airs_rad_gran_t * gran,
  double *nu,
  int nd,
  wave_t * wave) {

  double x0[3], x1[3];

  int ichan[AIRS_RAD_CHANNEL], id, track, xtrack;

  /* Get channel numbers... */
  for (id = 0; id < nd; id++) {
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
  for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
    geo2cart(0, gran->Longitude[0][xtrack], gran->Latitude[0][xtrack], x1);
    wave->x[xtrack] = DIST(x0, x1);
  }
  for (track = 0; track < AIRS_RAD_GEOTRACK; track++) {
    geo2cart(0, gran->Longitude[track][0], gran->Latitude[track][0], x1);
    wave->y[track] = DIST(x0, x1);
  }

  /* Set geolocation... */
  wave->time =
    gran->Time[AIRS_RAD_GEOTRACK / 2][AIRS_RAD_GEOXTRACK / 2] - 220838400;
  wave->z = 0;
  for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
    for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
      wave->lon[xtrack][track] = gran->Longitude[track][xtrack];
      wave->lat[xtrack][track] = gran->Latitude[track][xtrack];
    }

  /* Set brightness temperature... */
  for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
    for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
      wave->temp[xtrack][track] = 0;
      wave->bg[xtrack][track] = 0;
      wave->pt[xtrack][track] = 0;
      wave->var[xtrack][track] = 0;
      for (id = 0; id < nd; id++) {
	if ((gran->state[track][xtrack] != 0)
	    || (gran->ExcludedChans[ichan[id]] > 2)
	    || (gran->CalChanSummary[ichan[id]] & 8)
	    || (gran->CalChanSummary[ichan[id]] & (32 + 64))
	    || (gran->CalFlag[track][ichan[id]] & 16))
	  wave->temp[xtrack][track] = GSL_NAN;
	else
	  wave->temp[xtrack][track]
	    += brightness(gran->radiances[track][xtrack][ichan[id]] * 1e-3,
			  gran->nominal_freq[ichan[id]]) / nd;
      }
    }
}

/*****************************************************************************/

void ret2wave(
  ret_t * ret,
  wave_t * wave,
  int dataset,
  int ip) {

  double x0[3], x1[3];

  int ids, ix, iy;

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
  for (ids = 0; ids < ret->nds; ids++) {

    /* Get horizontal indices... */
    ix = ids % 90;
    iy = ids / 90;

    /* Get distances... */
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
  wave_t * wave,
  double dh) {

  double dh2, mu, help;

  int dx, dy, ix, ix2, iy, iy2, n;

  /* Check parameters... */
  if (dh <= 0)
    return;

  /* Compute squared radius... */
  dh2 = gsl_pow_2(dh);

  /* Get sampling distances... */
  dx =
    (int) (dh / fabs(wave->x[wave->nx - 1] - wave->x[0]) * (wave->nx - 1.0) +
	   1);
  dy =
    (int) (dh / fabs(wave->y[wave->ny - 1] - wave->y[0]) * (wave->ny - 1.0) +
	   1);

  /* Loop over data points... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++) {

      /* Init... */
      mu = help = 0;
      n = 0;

      /* Get data... */
      for (ix2 = GSL_MAX(ix - dx, 0); ix2 <= GSL_MIN(ix + dx, wave->nx - 1);
	   ix2++)
	for (iy2 = GSL_MAX(iy - dy, 0); iy2 <= GSL_MIN(iy + dy, wave->ny - 1);
	     iy2++)
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
  airs_l1_t * l1) {

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
  airs_l2_t * l2) {

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
  wave_t * wave) {

  FILE *out;

  int i, j;

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
	  "# $9  = perturbation [K]\n" "# $10 = variance [K^2]\n");

  /* Write data... */
  for (j = 0; j < wave->ny; j++) {
    fprintf(out, "\n");
    for (i = 0; i < wave->nx; i++)
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g\n",
	      wave->time, wave->z, wave->lon[i][j], wave->lat[i][j],
	      wave->x[i], wave->y[j], wave->temp[i][j], wave->bg[i][j],
	      wave->pt[i][j], wave->var[i][j]);
  }

  /* Close file... */
  fclose(out);
}
