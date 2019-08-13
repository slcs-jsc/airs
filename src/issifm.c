#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum model dimensions (ICON).
   #define NLON 1751
   #define NLAT 1201
   #define NZ 242
*/

/* Maximum model dimensions (IFS).
   #define NLON 1441
   #define NLAT 721
   #define NZ 138
*/

/* Maximum model dimensions (UM). */
#define NLON 3000
#define NLAT 910
#define NZ 185

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Interpolation of model data. */
void intpol(
  float ps[NLON][NLAT][NZ],
  float ts[NLON][NLAT][NZ],
  float zs[NLON][NLAT][NZ],
  double lons[NLON],
  double lats[NLAT],
  int nz,
  int nlon,
  int nlat,
  double z,
  double lon,
  double lat,
  double *p,
  double *t);

/*! Smoothing of model data. */
void smooth(
  float ps[NLON][NLAT][NZ],
  float ts[NLON][NLAT][NZ],
  float zs[NLON][NLAT][NZ],
  double lons[NLON],
  double lats[NLAT],
  int nz,
  int nlon,
  int nlat);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;

  static char kernel[LEN], pertname[LEN];

  static double z_um[182]
    = { 5.52994, 5.70933, 5.89159, 6.07672, 6.2647, 6.45555, 6.64927, 6.84585,
    7.04529, 7.2476, 7.45277, 7.6608, 7.8717, 8.08546, 8.30209, 8.52157,
    8.74393, 8.96914, 9.19722, 9.42817, 9.66198, 9.89865, 10.1382, 10.3806,
    10.6258, 10.874, 11.125, 11.3788, 11.6355, 11.8951, 12.1576, 12.4229,
    12.691, 12.9621, 13.236, 13.5128, 13.7924, 14.0749, 14.3602, 14.6485,
    14.9395, 15.2335, 15.5303, 15.83, 16.1325, 16.4379, 16.7462, 17.0573,
    17.3713, 17.6882, 18.0079, 18.3305, 18.6559, 18.9842, 19.3154, 19.6494,
    19.9864, 20.3261, 20.6688, 21.0143, 21.3626, 21.7138, 22.0679, 22.4249,
    22.7847, 23.1474, 23.5129, 23.8813, 24.2526, 24.6267, 25.0037, 25.3836,
    25.7663, 26.1519, 26.5403, 26.9317, 27.3258, 27.7229, 28.1228, 28.5256,
    28.9312, 29.3397, 29.7511, 30.1653, 30.5824, 31.0023, 31.4252, 31.8508,
    32.2794, 32.7108, 33.1451, 33.5822, 34.0222, 34.4651, 34.9108, 35.3594,
    35.8109, 36.2652, 36.7224, 37.1824, 37.6453, 38.1111, 38.5797, 39.0512,
    39.5256, 40.0028, 40.4829, 40.9659, 41.4517, 41.9404, 42.4319, 42.9264,
    43.4236, 43.9238, 44.4268, 44.9327, 45.4414, 45.953, 46.4674, 46.9848,
    47.505, 48.028, 48.5539, 49.0827, 49.6143, 50.1488, 50.6862, 51.2265,
    51.7696, 52.3155, 52.8643, 53.416, 53.9706, 54.528, 55.0883, 55.6514,
    56.2174, 56.7863, 57.358, 57.9326, 58.5101, 59.0904, 59.6736, 60.2597,
    60.8486, 61.4404, 62.035, 62.6325, 63.2329, 63.8361, 64.4423, 65.0512,
    65.663, 66.2777, 66.8953, 67.5157, 68.139, 68.7651, 69.3942, 70.026,
    70.6608, 71.2984, 71.9388, 72.5822, 73.2284, 73.8774, 74.5294, 75.1841,
    75.8418, 76.5023, 77.1657, 77.8319, 78.501, 79.173, 79.8478, 80.5255,
    81.2061, 81.8895, 82.5758, 83.2649, 83.957, 84.6518
  };

  static double lon[NLON], lat[NLAT], xo[3], xs[3], xm[3], var_dh = 100.,
    f, t_ovp, hyam[NZ], hybm[NZ], kz[NSHAPE], kw[NSHAPE], w, wsum;

  static float *help, ps[NLON][NLAT], p[NLON][NLAT][NZ], t[NLON][NLAT][NZ],
    z[NLON][NLAT][NZ];

  static int init, id, itrack, ixtrack, ncid, dimid, varid, slant,
    ilon, ilat, iz, nlon, nlat, nz, nz2, ip, track0, track1, nk, okay;

  static size_t rs;

  atm_t *atm;

  obs_t *obs;

  pert_t *pert;

  wave_t *wave;

  /* ------------------------------------------------------------
     Get control parameters...
     ------------------------------------------------------------ */

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <model> <model.nc> <pert.nc>"
	   " <wave_airs.tab> <wave_model.tab>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  scan_ctl(argc, argv, "KERNEL", -1, "-", kernel);
  slant = (int) scan_ctl(argc, argv, "SLANT", -1, "1", NULL);
  t_ovp = scan_ctl(argc, argv, "T_OVP", -1, "", NULL);

  /* Set control parameters... */
  ctl.write_bbt = 1;

  /* ------------------------------------------------------------
     Read model data...
     ------------------------------------------------------------ */

  /* Allocate... */
  ALLOC(help, float,
	NLON * NLAT * NZ);

  /* Read ICON data... */
  if (strcasecmp(argv[2], "icon") == 0) {

    /* Open file... */
    printf("Read ICON data: %s\n", argv[3]);
    NC(nc_open(argv[3], NC_NOWRITE, &ncid));

    /* Get dimensions... */
    NC(nc_inq_dimid(ncid, "height", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nz = (int) rs;
    if (nz > NZ)
      ERRMSG("Too many altitudes!");

    NC(nc_inq_dimid(ncid, "lat", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nlat = (int) rs;
    if (nlat > NLAT)
      ERRMSG("Too many latitudes!");

    NC(nc_inq_dimid(ncid, "lon", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nlon = (int) rs;
    if (nlon > NLON)
      ERRMSG("Too many longitudes!");

    /* Read latitudes... */
    NC(nc_inq_varid(ncid, "lat", &varid));
    NC(nc_get_var_double(ncid, varid, lat));

    /* Read longitudes... */
    NC(nc_inq_varid(ncid, "lon", &varid));
    NC(nc_get_var_double(ncid, varid, lon));

    /* Read temperature... */
    NC(nc_inq_varid(ncid, "temp", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  t[ilon][ilat][iz] = help[(iz * nlat + ilat) * nlon + ilon];

    /* Read height... */
    NC(nc_inq_varid(ncid, "z_mc", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  z[ilon][ilat][iz] =
	    (float) (help[(iz * nlat + ilat) * nlon + ilon] / 1e3);

    /* Read pressure... */
    NC(nc_inq_varid(ncid, "pres", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  p[ilon][ilat][iz] =
	    (float) (help[(iz * nlat + ilat) * nlon + ilon] / 1e2);

    /* Close file... */
    NC(nc_close(ncid));
  }

  /* Read IFS data... */
  else if (strcasecmp(argv[2], "ifs") == 0) {

    /* Open file... */
    printf("Read IFS data: %s\n", argv[3]);
    NC(nc_open(argv[3], NC_NOWRITE, &ncid));

    /* Get dimensions... */
    NC(nc_inq_dimid(ncid, "lev_2", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nz = (int) rs;
    if (nz > NZ)
      ERRMSG("Too many altitudes!");

    NC(nc_inq_dimid(ncid, "lat", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nlat = (int) rs;
    if (nlat > NLAT)
      ERRMSG("Too many latitudes!");

    NC(nc_inq_dimid(ncid, "lon", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nlon = (int) rs;
    if (nlon > NLON)
      ERRMSG("Too many longitudes!");

    /* Read latitudes... */
    NC(nc_inq_varid(ncid, "lat", &varid));
    NC(nc_get_var_double(ncid, varid, lat));

    /* Read longitudes... */
    NC(nc_inq_varid(ncid, "lon", &varid));
    NC(nc_get_var_double(ncid, varid, lon));

    /* Read temperature... */
    NC(nc_inq_varid(ncid, "t", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  t[ilon][ilat][iz] = help[(iz * nlat + ilat) * nlon + ilon];

    /* Read height... */
    NC(nc_inq_varid(ncid, "gh", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  z[ilon][ilat][iz] =
	    (float) (help[(iz * nlat + ilat) * nlon + ilon] / 1e3);

    /* Read surface pressure... */
    NC(nc_inq_varid(ncid, "lnsp", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	ps[ilon][ilat] = (float) exp(help[ilat * nlon + ilon]);

    /* Read grid coefficients... */
    NC(nc_inq_varid(ncid, "hyam", &varid));
    NC(nc_get_var_double(ncid, varid, hyam));
    NC(nc_inq_varid(ncid, "hybm", &varid));
    NC(nc_get_var_double(ncid, varid, hybm));

    /* Calculate pressure... */
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  p[ilon][ilat][iz]
	    = (float) ((hyam[iz] + hybm[iz] * ps[ilon][ilat]) / 100.);

    /* Close file... */
    NC(nc_close(ncid));
  }

  /* Read UM data... */
  else if (strcasecmp(argv[2], "um") == 0) {

    /* Open file... */
    printf("Read UM data: %s\n", argv[3]);
    NC(nc_open(argv[3], NC_NOWRITE, &ncid));

    /* Get dimensions... */
    NC(nc_inq_dimid(ncid, "RHO_TOP_eta_rho", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nz = (int) rs;
    if (nz > NZ)
      ERRMSG("Too many altitudes!");
    NC(nc_inq_dimid(ncid, "RHO_BOT_eta_rho", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nz2 = (int) rs;
    if (nz2 > NZ)
      ERRMSG("Too many altitudes!");

    NC(nc_inq_dimid(ncid, "latitude", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nlat = (int) rs;
    if (nlat > NLAT)
      ERRMSG("Too many latitudes!");

    NC(nc_inq_dimid(ncid, "longitude", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    nlon = (int) rs;
    if (nlon > NLON)
      ERRMSG("Too many longitudes!");

    /* Read latitudes... */
    NC(nc_inq_varid(ncid, "latitude", &varid));
    NC(nc_get_var_double(ncid, varid, lat));

    /* Read longitudes... */
    NC(nc_inq_varid(ncid, "longitude", &varid));
    NC(nc_get_var_double(ncid, varid, lon));

    /* Read temperature... */
    NC(nc_inq_varid(ncid, "STASH_m01s30i004", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  t[ilon][ilat][iz] = help[(iz * nlat + ilat) * nlon + ilon];

    /* Read pressure... */
    NC(nc_inq_varid(ncid, "STASH_m01s00i407", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  p[ilon][ilat][iz] = 0.01f * help[(iz * nlat + ilat) * nlon + ilon];

    /* Set fixed height levels... */
    if (nz != 182)
      ERRMSG("Wrong number of height levels!");
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  z[ilon][ilat][iz] = (float) z_um[iz];

    /* Read low-level heights... */
    NC(nc_inq_varid(ncid, "STASH_m01s15i102", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz2; iz++)
	  z[ilon][ilat][iz] = (float) (help[iz] / 1e3);

    /* Check data... */
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  if (t[ilon][ilat][iz] <= 100 || t[ilon][ilat][iz] >= 400) {
	    p[ilon][ilat][iz] = GSL_NAN;
	    t[ilon][ilat][iz] = GSL_NAN;
	  }

    /* Close file... */
    NC(nc_close(ncid));
  }

  else
    ERRMSG("Model type not supported!");

  /* Free... */
  free(help);

  /* Smoothing of model data... */
  smooth(p, t, z, lon, lat, nz, nlon, nlat);

  /* ------------------------------------------------------------
     Read AIRS perturbation data...
     ------------------------------------------------------------ */

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(obs, obs_t, 1);
  ALLOC(pert, pert_t, 1);
  ALLOC(wave, wave_t, 1);

  /* Read perturbation data... */
  read_pert(argv[4], pertname, pert);

  /* Find track range... */
  for (itrack = 0; itrack < pert->ntrack; itrack++) {
    if (pert->time[itrack][44] < t_ovp - 720 || itrack == 0)
      track0 = itrack;
    track1 = itrack;
    if (pert->time[itrack][44] > t_ovp + 720)
      break;
  }

  /* Convert to wave analysis struct... */
  pert2wave(pert, wave, track0, track1, 0, pert->nxtrack - 1);

  /* Estimate background... */
  background_poly(wave, 5, 0);

  /* Compute variance... */
  variance(wave, var_dh);

  /* Write observation wave struct... */
  write_wave(argv[5], wave);

  /* ------------------------------------------------------------
     Run forward model...
     ------------------------------------------------------------ */

  /* Loop over AIRS geolocations... */
  for (itrack = track0; itrack <= track1; itrack++)
    for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {

      /* Write info... */
      if (ixtrack == 0)
	printf("Compute track %d / %d ...\n", itrack - track0 + 1,
	       track1 - track0 + 1);

      /* Set observation data... */
      obs->nr = 1;
      obs->obsz[0] = 705;
      obs->obslon[0] = pert->lon[itrack][44];
      obs->obslat[0] = pert->lat[itrack][44];

      /* Get Cartesian coordinates... */
      geo2cart(obs->obsz[0], obs->obslon[0], obs->obslat[0], xo);
      geo2cart(0, pert->lon[itrack][ixtrack], pert->lat[itrack][ixtrack], xs);

      /* Set profile for atmospheric data... */
      if (slant) {
	atm->np = 0;
	for (f = 0.0; f <= 1.0; f += 0.0002) {
	  xm[0] = f * xo[0] + (1 - f) * xs[0];
	  xm[1] = f * xo[1] + (1 - f) * xs[1];
	  xm[2] = f * xo[2] + (1 - f) * xs[2];
	  cart2geo(xm, &atm->z[atm->np], &atm->lon[atm->np],
		   &atm->lat[atm->np]);
	  atm->time[atm->np] = pert->time[itrack][ixtrack];
	  if (atm->z[atm->np] < 10)
	    continue;
	  else if (atm->z[atm->np] > 90)
	    break;
	  else if ((++atm->np) >= NP)
	    ERRMSG("Too many altitudes!");
	}
      } else {
	atm->np = 0;
	for (f = 10.0; f <= 90.0; f += 0.2) {
	  atm->time[atm->np] = pert->time[itrack][ixtrack];
	  atm->z[atm->np] = f;
	  atm->lon[atm->np] = pert->lon[itrack][ixtrack];
	  atm->lat[atm->np] = pert->lat[itrack][ixtrack];
	  if ((++atm->np) >= NP)
	    ERRMSG("Too many altitudes!");
	}
      }

      /* Initialize with climatological data... */
      climatology(&ctl, atm);

      /* Interpolate model data... */
      for (ip = 0; ip < atm->np; ip++)
	intpol(p, t, z, lon, lat, nz, nlon, nlat, atm->z[ip],
	       atm->lon[ip], atm->lat[ip], &atm->p[ip], &atm->t[ip]);

      /* Check profile... */
      okay = 1;
      for (ip = 0; ip < atm->np; ip++)
	if (!gsl_finite(atm->p[ip]) || !gsl_finite(atm->t[ip]))
	  okay = 0;
      if (!okay)
	pert->bt[itrack][ixtrack] = GSL_NAN;
      else {

	/* Use kernel function... */
	if (kernel[0] != '-') {

	  /* Read kernel function... */
	  if (!init) {
	    init = 1;
	    read_shape(kernel, kz, kw, &nk);
	    if (kz[0] > kz[1])
	      ERRMSG("Kernel function must be ascending!");
	  }

	  /* Calculate mean temperature... */
	  pert->bt[itrack][ixtrack] = wsum = 0;
	  for (ip = 0; ip < atm->np; ip++)
	    if (atm->z[ip] >= kz[0] && atm->z[ip] <= kz[nk - 1]) {
	      iz = locate_irr(kz, nk, atm->z[ip]);
	      w = LIN(kz[iz], kw[iz], kz[iz + 1], kw[iz + 1], atm->z[ip]);
	      pert->bt[itrack][ixtrack] += w * atm->t[ip];
	      wsum += w;
	    }
	  pert->bt[itrack][ixtrack] /= wsum;
	}

	/* Use radiative transfer model... */
	else {

	  /* Run forward model... */
	  formod(&ctl, atm, obs);

	  /* Get mean brightness temperature... */
	  pert->bt[itrack][ixtrack] = 0;
	  for (id = 0; id < ctl.nd; id++)
	    pert->bt[itrack][ixtrack] += obs->rad[id][0] / ctl.nd;
	}
      }
    }

  /* ------------------------------------------------------------
     Write model perturbations...
     ------------------------------------------------------------ */

  /* Convert to wave analysis struct... */
  pert2wave(pert, wave, track0, track1, 0, pert->nxtrack - 1);

  /* Estimate background... */
  background_poly(wave, 5, 0);

  /* Compute variance... */
  variance(wave, var_dh);

  /* Write observation wave struct... */
  write_wave(argv[6], wave);

  /* Free... */
  free(atm);
  free(obs);
  free(pert);
  free(wave);

  return EXIT_SUCCESS;
}

/************************************************************************/

void intpol(
  float ps[NLON][NLAT][NZ],
  float ts[NLON][NLAT][NZ],
  float zs[NLON][NLAT][NZ],
  double lons[NLON],
  double lats[NLAT],
  int nz,
  int nlon,
  int nlat,
  double z,
  double lon,
  double lat,
  double *p,
  double *t) {

  double p00, p01, p10, p11, t00, t01, t10, t11, zd[NZ];

  int iz, ilon, ilat;

  /* Adjust longitude... */
  if (lons[nlon - 1] > 180)
    if (lon < 0)
      lon += 360;

  /* Get indices... */
  ilon = locate_reg(lons, nlon, lon);
  ilat = locate_reg(lats, nlat, lat);

  /* Check vertical range... */
  if (z > GSL_MAX(zs[ilon][ilat][0], zs[ilon][ilat][nz - 1])
      || z < GSL_MIN(zs[ilon][ilat][0], zs[ilon][ilat][nz - 1])
      || z > GSL_MAX(zs[ilon][ilat + 1][0], zs[ilon][ilat + 1][nz - 1])
      || z < GSL_MIN(zs[ilon][ilat + 1][0], zs[ilon][ilat + 1][nz - 1])
      || z > GSL_MAX(zs[ilon + 1][ilat][0], zs[ilon + 1][ilat][nz - 1])
      || z < GSL_MIN(zs[ilon + 1][ilat][0], zs[ilon + 1][ilat][nz - 1])
      || z > GSL_MAX(zs[ilon + 1][ilat + 1][0],
		     zs[ilon + 1][ilat + 1][nz - 1])
      || z < GSL_MIN(zs[ilon + 1][ilat + 1][0],
		     zs[ilon + 1][ilat + 1][nz - 1]))
    return;

  /* Interpolate vertically... */
  for (iz = 0; iz < nz; iz++)
    zd[iz] = zs[ilon][ilat][iz];
  iz = locate_irr(zd, nz, z);
  p00 = LIN(zs[ilon][ilat][iz], ps[ilon][ilat][iz],
	    zs[ilon][ilat][iz + 1], ps[ilon][ilat][iz + 1], z);
  t00 = LIN(zs[ilon][ilat][iz], ts[ilon][ilat][iz],
	    zs[ilon][ilat][iz + 1], ts[ilon][ilat][iz + 1], z);

  for (iz = 0; iz < nz; iz++)
    zd[iz] = zs[ilon][ilat + 1][iz];
  iz = locate_irr(zd, nz, z);
  p01 = LIN(zs[ilon][ilat + 1][iz], ps[ilon][ilat + 1][iz],
	    zs[ilon][ilat + 1][iz + 1], ps[ilon][ilat + 1][iz + 1], z);
  t01 = LIN(zs[ilon][ilat + 1][iz], ts[ilon][ilat + 1][iz],
	    zs[ilon][ilat + 1][iz + 1], ts[ilon][ilat + 1][iz + 1], z);

  for (iz = 0; iz < nz; iz++)
    zd[iz] = zs[ilon + 1][ilat][iz];
  iz = locate_irr(zd, nz, z);
  p10 = LIN(zs[ilon + 1][ilat][iz], ps[ilon + 1][ilat][iz],
	    zs[ilon + 1][ilat][iz + 1], ps[ilon + 1][ilat][iz + 1], z);
  t10 = LIN(zs[ilon + 1][ilat][iz], ts[ilon + 1][ilat][iz],
	    zs[ilon + 1][ilat][iz + 1], ts[ilon + 1][ilat][iz + 1], z);

  for (iz = 0; iz < nz; iz++)
    zd[iz] = zs[ilon + 1][ilat + 1][iz];
  iz = locate_irr(zd, nz, z);
  p11 = LIN(zs[ilon + 1][ilat + 1][iz], ps[ilon + 1][ilat + 1][iz],
	    zs[ilon + 1][ilat + 1][iz + 1], ps[ilon + 1][ilat + 1][iz + 1],
	    z);
  t11 = LIN(zs[ilon + 1][ilat + 1][iz], ts[ilon + 1][ilat + 1][iz],
	    zs[ilon + 1][ilat + 1][iz + 1], ts[ilon + 1][ilat + 1][iz + 1],
	    z);

  /* Interpolate horizontally... */
  p00 = LIN(lons[ilon], p00, lons[ilon + 1], p10, lon);
  p11 = LIN(lons[ilon], p01, lons[ilon + 1], p11, lon);
  *p = LIN(lats[ilat], p00, lats[ilat + 1], p11, lat);

  t00 = LIN(lons[ilon], t00, lons[ilon + 1], t10, lon);
  t11 = LIN(lons[ilon], t01, lons[ilon + 1], t11, lon);
  *t = LIN(lats[ilat], t00, lats[ilat + 1], t11, lat);
}

/************************************************************************/

void smooth(
  float ps[NLON][NLAT][NZ],
  float ts[NLON][NLAT][NZ],
  float zs[NLON][NLAT][NZ],
  double lons[NLON],
  double lats[NLAT],
  int nz,
  int nlon,
  int nlat) {

  static float hp[NLON][NLAT], ht[NLON][NLAT], hz[NLON][NLAT], w, wsum;

  static double dx, dy, wx[10], wy[10];

  int iz, ilon, ilon2, ilat, ilat2, dlon = 3, dlat = 3;

  /* Set weights... */
  dy = RE * M_PI / 180. * fabs(lats[1] - lats[0]);
  for (ilat = 0; ilat <= dlat; ilat++)
    wy[ilat] = exp(-0.5 * POW2(ilat * dy * 2.35482 / 20.));

  /* Loop over height levels... */
  for (iz = 0; iz < nz; iz++) {

    /* Write info... */
    printf("Smoothing level %d / %d ...\n", iz + 1, nz);

    /* Copy data... */
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++) {
	hp[ilon][ilat] = ps[ilon][ilat][iz];
	ht[ilon][ilat] = ts[ilon][ilat][iz];
	hz[ilon][ilat] = zs[ilon][ilat][iz];
      }

    /* Loop over latitudes... */
    for (ilat = 0; ilat < nlat; ilat++) {

      /* Set weights... */
      dx = RE * M_PI / 180. * cos(lats[ilat] * M_PI / 180.) *
	fabs(lons[1] - lons[0]);
      for (ilon = 0; ilon <= dlon; ilon++)
	wx[ilon] = exp(-0.5 * POW2(ilon * dx * 2.35482 / 20.));

      /* Loop over longitudes... */
      for (ilon = 0; ilon < nlon; ilon++) {
	wsum = 0;
	ps[ilon][ilat][iz] = 0;
	ts[ilon][ilat][iz] = 0;
	zs[ilon][ilat][iz] = 0;
	for (ilon2 = GSL_MAX(ilon - dlon, 0);
	     ilon2 <= GSL_MIN(ilon + dlon, nlon - 1); ilon2++)
	  for (ilat2 = GSL_MAX(ilat - dlat, 0);
	       ilat2 <= GSL_MIN(ilat + dlat, nlat - 1); ilat2++) {
	    w = (float) (wx[abs(ilon2 - ilon)] * wy[abs(ilat2 - ilat)]);
	    ps[ilon][ilat][iz] += w * hp[ilon2][ilat2];
	    ts[ilon][ilat][iz] += w * ht[ilon2][ilat2];
	    zs[ilon][ilat][iz] += w * hz[ilon2][ilat2];
	    wsum += w;
	  }
	ps[ilon][ilat][iz] /= wsum;
	ts[ilon][ilat][iz] /= wsum;
	zs[ilon][ilat][iz] /= wsum;
      }
    }
  }
}
