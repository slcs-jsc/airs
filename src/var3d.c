#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static atm_t atm;
  static obs_t obs;

  gsl_multifit_linear_workspace *work;
  gsl_matrix *cov, *k, *X;
  gsl_vector *c, *xvec, *yvec;

  static double alpha, alphamax, amp, ampmax, bg[L1_NXTRACK], ca, chisq,
    dalpha, ddx, dx, jac[L1_NXTRACK][NP], ly, lz, mu, phi, rad[L1_NXTRACK],
    radius, sa, t30, var, vmean, vmin, vmax, x, y[L1_NXTRACK];

  static int detrend, dim = 5, i, i2, id, ip, n, nmu, nphi, ndx;

  /* Check arguments... */
  if (argc < 13)
    ERRMSG("Give parameters: <ctl> <atm> <T0_30km> <exp/lin> <radius> "
	   "<obsz> <alphamax> <n> <dalpha> <dx> <ddx> <detrend>");
  t30 = atof(argv[3]);
  radius = atof(argv[5]);
  obs.obsz[0] = atof(argv[6]);
  alphamax = atof(argv[7]);
  n = atoi(argv[8]);
  if (n > L1_NXTRACK)
    ERRMSG("Too many tracks!");
  dalpha = atof(argv[9]);
  dx = atof(argv[10]);
  ddx = atof(argv[11]);
  detrend = atoi(argv[12]);

  /* Initialize... */
  c = gsl_vector_alloc((size_t) dim);
  cov = gsl_matrix_alloc((size_t) dim, (size_t) dim);
  work = gsl_multifit_linear_alloc((size_t) n, (size_t) dim);
  X = gsl_matrix_alloc((size_t) n, (size_t) dim);
  xvec = gsl_vector_alloc((size_t) n);
  yvec = gsl_vector_alloc((size_t) n);

  /* Read forward model control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read atmospheric data... */
  read_atm(NULL, argv[2], &ctl, &atm);

  /* ------------------------------------------------------------
     Compute mean radiance and kernel functions...
     ------------------------------------------------------------ */

  /* Loop over scans... */
  for (i = 0; i < n; i++) {

    /* Set observation geometry... */
    obs.nr = 1;
    alpha =
      -alphamax + 2. * alphamax * i / (n - 1.) + (i % 2 ==
						  0 ? 1.0 : -1.0) * dalpha;
    sa = sin(alpha * M_PI / 180.);
    ca = cos(alpha * M_PI / 180.);
    obs.vplat[0] = 180. / M_PI
      * asin(sa / RE * ((RE + obs.obsz[0]) * ca
			- sqrt(gsl_pow_2(RE) -
			       gsl_pow_2((RE + obs.obsz[0]) * sa))));
    y[i] = obs.vplat[0] / 180 * M_PI * RE;

    /* Run forward model... */
    formod(&ctl, &atm, &obs);
    bg[i] = 0;
    for (id = 0; id < ctl.nd; id++)
      bg[i] += obs.rad[id][0] / ctl.nd;

    /* Compute kernel matrix... */
    ctl.rett_zmin = -10000;
    ctl.rett_zmax = 10000;
    k = gsl_matrix_alloc((size_t) ctl.nd, (size_t) atm.np);
    kernel(&ctl, &atm, &obs, k);
    for (ip = 0; ip < atm.np; ip++) {
      jac[i][ip] = 0;
      for (id = 0; id < ctl.nd; id++)
	jac[i][ip] += gsl_matrix_get(k, (size_t) id, (size_t) ip) / ctl.nd;
    }
    gsl_matrix_free(k);
  }

  /* ------------------------------------------------------------
     Get variance filter characteristics...
     ------------------------------------------------------------ */

  /* Loop over wavelengths... */
  for (lz = 10; lz <= 50; lz += 0.5)
    for (ly = 50; ly <= 1500; ly += 10) {

      /* Initialize... */
      vmean = 0;
      vmin = 1e10;
      vmax = -1e10;
      nphi = 0;

      /* Loop over phases... */
      for (phi = 0; phi < 2 * M_PI; phi += M_PI / 24) {

	/* Initialize... */
	nmu = 0;
	mu = var = 0;

	/* Loop over swaths... */
	for (x = -radius; x <= radius;
	     x += dx + ((ndx++) % 2 == 0 ? 1.0 : -1.0) * ddx) {

	  /* Compute radiances for perturbed profile... */
	  for (i = 0; i < n; i++) {
	    rad[i] = bg[i];
	    for (ip = 0; ip < atm.np; ip++) {
	      amp = t30;
	      if (argv[4][0] == 'e' || argv[4][0] == 'E') {

		/* Saturation amplitude (Preusse et al., 2008),
		   Tmax = lz / (2*pi) * Tbg / g * N^2... */
		ampmax = lz * 1e3 / (2 * M_PI) * 250 / 9.81 * gsl_pow_2(0.02);

		/* Get wave amplitude... */
		amp *= exp((atm.z[ip] - 30.) / 14.);
		amp = (amp > ampmax) ? ampmax : amp;
	      }
	      rad[i] += jac[i][ip] * amp
		* sin(2 * M_PI / ly * y[i] + 2 * M_PI / lz * atm.z[ip] + phi);
	    }
	  }

	  /* Detrending... */
	  if (detrend) {
	    for (i = 0; i < n; i++) {
	      gsl_vector_set(xvec, (size_t) i, y[i]);
	      gsl_vector_set(yvec, (size_t) i, rad[i]);
	      for (i2 = 0; i2 < dim; i2++)
		gsl_matrix_set(X, (size_t) i, (size_t) i2,
			       pow(gsl_vector_get(xvec, (size_t) i),
				   1. * i2));
	    }
	    gsl_multifit_linear(X, yvec, c, cov, &chisq, work);
	    for (i = 0; i < n; i++)
	      rad[i] -= gsl_poly_eval(c->data, (int) dim,
				      gsl_vector_get(xvec, (size_t) i));
	  }

	  /* Compute variance... */
	  for (i = 0; i < n; i++)
	    if (gsl_pow_2(x) + gsl_pow_2(y[i]) <= gsl_pow_2(radius)) {
	      mu += rad[i];
	      var += gsl_pow_2(rad[i]);
	      nmu++;
	    }
	}

	/* Compute variance... */
	mu /= nmu;
	var = var / nmu - mu * mu;
	vmean += var;
	vmax = GSL_MAX(vmax, var);
	vmin = GSL_MIN(vmin, var);
	nphi++;
      }

      /* Write output... */
      printf("obsfilt: %g %g %g %g %g\n", ly, lz, vmean / nphi, vmax, vmin);
    }

  return EXIT_SUCCESS;
}
