#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of levels. */
#define NZ 1000

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Compute buoyancy frequency. */
double buoyancy(
  double z0,
  double p0,
  double t0,
  double z1,
  double p1,
  double t1);

/* Compute scale height. */
double scale_height(
  double t);

/* Convert temperature to potential temperature. */
double temp2theta(
  double p,
  double t);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  FILE *in;

  static double f0, k, omin, z[NZ], u[NZ], urel[NZ], v[NZ], bf[NZ], bf2[NZ],
    H[NZ], frel[NZ], osign[NZ], f1[NZ], f2[NZ], delta[NZ], a2[NZ], m[NZ],
    dxdz[NZ], cgz[NZ], dz, path[NZ], tim[NZ], costh, p[NZ], t[NZ], z0, w,
    wsum, dzw = 5 * 1e3, fgb, m0, alpha, lat;

  static int iz, iz2, izcrit, izrefl, nz;

  /* Check arguments... */
  if (argc != 8)
    ERRMSG("Give parameters: <atm.tab> <z_launch> <mode> "
	   "<t_gb | lz_launch> <lx> <lat> <direct>");

  /* Get launch level... */
  z0 = atof(argv[2]);
  lat = atof(argv[6]);
  alpha = atof(argv[7]);

  /* Read atmosphere above launch level... */
  if (!(in = fopen(argv[1], "r")))
    ERRMSG("Cannot open atmospheric data file!");
  while (fscanf
	 (in, "%lg %lg %lg %lg %lg", &z[nz], &p[nz], &t[nz], &u[nz], &v[nz])
	 == 5)
    if (z[nz] >= z0) {
      u[nz] =
	cos(alpha * M_PI / 180.) * u[nz] + sin(alpha * M_PI / 180.) * v[nz];
      if ((++nz) > NZ)
	ERRMSG("Too many altitude levels!");
    }
  fclose(in);

  /* Compute scale height and buoyancy frequency... */
  for (iz = 0; iz < nz; iz++) {
    if (iz < nz - 1)
      bf[iz] = buoyancy(z[iz], p[iz], t[iz], z[iz + 1], p[iz + 1], t[iz + 1]);
    else
      bf[iz] = bf[iz - 1];
    H[iz] = scale_height(t[iz]) * 1e3;
    z[iz] *= 1e3;
  }

  /* Smooth N profile... */
  for (iz = 0; iz < nz; iz++) {
    bf2[iz] = wsum = 0;
    for (iz2 = 0; iz2 < nz; iz2++) {
      if (!gsl_finite(bf[iz2]) ||
	  !gsl_finite(bf[GSL_MAX(iz2 - 1, 0)]) ||
	  !gsl_finite(bf[GSL_MIN(iz2 + 1, nz - 1)]))
	continue;
      w =
	(fabs(z[iz] - z[iz2]) < dzw) ? 1.0 - fabs(z[iz] - z[iz2]) / dzw : 0.0;
      bf2[iz] += w * bf[iz2];
      wsum += w;
    }
    bf2[iz] /= wsum;
  }
  for (iz = 0; iz < nz; iz++)
    bf[iz] = bf2[iz];

  /* Get horizontal wavenumber... */
  k = 2 * M_PI / (atof(argv[5]) * 1e3);

  /* Get minimum gravity wave frequency (Coriolis parameter)... */
  omin = 2 * 2 * M_PI / 86400. * sin(lat / 180. * M_PI);

  /* Get initial frequencies... */
  if (argv[3][0] == 't') {

    /* Get ground-based frequency... */
    fgb = 2 * M_PI / (atof(argv[4]) * 60.);

    /* Get intrinsic frequency at launch level... */
    f0 = fgb - k * u[0];

  } else if (argv[3][0] == 'l') {

    /* Get vertical wavenumber... */
    m0 = 2 * M_PI / (atof(argv[4]) * 1e3);

    /* Get intrinsic frequency at launch level... */
    f0 =
      sqrt((bf[0] * bf[0] * k * k +
	    omin * omin * (m0 * m0 + 0.25 / (H[0] * H[0])))
	   / (m0 * m0 + k * k + 0.25 / (H[0] * H[0])));

    /* Get ground-based frequency... */
    fgb = f0 + k * u[0];

  } else
    ERRMSG("Set <mode> to 't_gb' or 'lz_launch'!");

  /* Loop over layers... */
  for (iz = 0; iz < nz; iz++) {
    urel[iz] = u[iz] - u[0];
    frel[iz] = f0 - k * urel[iz];
    osign[iz] = frel[iz] / fabs(frel[iz]);
    f1[iz] = (bf[iz] * bf[iz] - frel[iz] * frel[iz]) / frel[iz];
    f2[iz] = (frel[iz] * frel[iz] - omin * omin) / frel[iz];
    delta[iz] = k * k * (1 + f1[iz] / f2[iz]);
    a2[iz] = 1. / 4. / (H[iz] * H[iz]);
    m[iz] = (-osign[iz]) * k * sqrt((f1[iz] / f2[iz]) - (a2[iz] / (k * k)));
    dxdz[iz] = (u[iz] * delta[iz] + k * f1[iz]) / (-1 * m[iz] * f2[iz]);
    dz = z[1] - z[0];
    cgz[iz] = f2[iz] * (-1. * m[iz]) / (k * k + m[iz] * m[iz] + a2[iz]);
  }

  /* Integrate via trapezoidal rule... */
  for (iz = 1; iz < nz; iz++) {
    path[iz] = path[iz - 1] + dz * .5 * (dxdz[iz - 1] + dxdz[iz]);
    tim[iz] = tim[iz - 1] + dz * 2. / (cgz[iz - 1] + cgz[iz]);
  }

  /* Find critical level... */
  for (izcrit = 0; izcrit < nz; izcrit++)
    if (f0 / fabs(f0) * frel[izcrit] / fabs(omin) <= 1)
      break;

  /* Find trapping/reflection level... */
  for (izrefl = 0; izrefl < nz; izrefl++) {
    costh = fabs(f0 - k * urel[izrefl])
      / sqrt(bf[izrefl] * bf[izrefl]
	     * (1 -
		(1 -
		 (omin / bf[izrefl]) * (omin / bf[izrefl])) / (k * k /
							       a2[izrefl] +
							       1)));
    if (costh >= 1.0)
      break;
  }

  /* Filter data... */
  for (iz = 0; iz < nz; iz++)
    if (iz >= izcrit || iz >= izrefl)
      path[iz] = tim[iz] = m[iz] = frel[iz] = cgz[iz] = sqrt(-1.0);

  /* Write output... */
  printf("# $1  = latitude [deg]\n"
	 "# $2  = altitude [km]\n"
	 "# $3  = pressure [hPa]\n"
	 "# $4  = temperature [K]\n"
	 "# $5  = potential temperature [K]\n"
	 "# $6  = wind speed [m/s]\n"
	 "# $7  = buoyancy frequency [1/s]\n"
	 "# $8  = scale height [km]\n"
	 "# $9  = horizontal distance [km]\n"
	 "# $10 = propagation time [min]\n"
	 "# $11 = vertical wavelength [km]\n"
	 "# $12 = wave period [min]\n"
	 "# $13 = vertical group velocity [m/s]\n\n");
  for (iz = 0; iz < nz; iz++)
    printf("%g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   lat, z[iz] / 1e3, p[iz], t[iz], temp2theta(p[iz], t[iz]), u[iz],
	   bf[iz], H[iz] / 1e3, path[iz] / 1e3, tim[iz] / 60,
	   fabs(2 * M_PI / m[iz] / 1e3), 2. * M_PI / frel[iz] / 60., cgz[iz]);
  printf("\n# z_crit= %g km\n# z_refl= %g km\n",
	 z[izcrit - 1] / 1e3, z[izrefl - 1] / 1e3);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

double buoyancy(
  double z0,
  double p0,
  double t0,
  double z1,
  double p1,
  double t1) {

  double theta0, theta1;

  /* Get potential temperature... */
  theta0 = temp2theta(p0, t0);
  theta1 = temp2theta(p1, t1);

  /* Get buoyancy frequency... */
  return sqrt(G0 / (0.5 * (theta0 + theta1)) * (theta1 - theta0) /
	      ((z1 - z0) * 1e3));
}

/*****************************************************************************/

double scale_height(
  double t) {

  return 29.26 * t / 1e3;
}

/*****************************************************************************/

double temp2theta(
  double p,
  double t) {

  return t * pow(P0 / p, 0.286);
}
