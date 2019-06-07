#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum model dimensions. */
#define NLON 1441
#define NLAT 721
#define NZ 138
#if 0
#define NLON 1751
#define NLAT 1201
#define NZ 242
#endif

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

  static double lon[NLON], lat[NLAT], xo[3], xs[3], xm[3], var_dh = 100.,
    f, t_ovp, hyam[NZ], hybm[NZ], kz[NSHAPE], kw[NSHAPE], w, wsum;

  static float *help, ps[NLON][NLAT], p[NLON][NLAT][NZ], t[NLON][NLAT][NZ],
    z[NLON][NLAT][NZ];

  static int init, id, itrack, ixtrack, ncid, dimid, varid, slant,
    ilon, ilat, iz, nlon, nlat, nz, ip, track0, track1, nk;

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

    /* Read geopotential heights... */
    NC(nc_inq_varid(ncid, "z_mc", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  z[ilon][ilat][iz] =
	    (float) (help[(iz * nlat + ilat) * nlon + ilon] / 1e3);

    /* Calculate pressure... */
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++)
	for (iz = 0; iz < nz; iz++)
	  p[ilon][ilat][iz]
	    = (float) (1013.25 * exp(-z[ilon][ilat][iz] / 7.0));

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

    /* Read geopotential heights... */
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
	intpol(p, t, z, lon, lat, nz, nlon, nlat,
	       atm->z[ip], atm->lon[ip], atm->lat[ip], &atm->p[ip],
	       &atm->t[ip]);

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
  if (z > zs[ilon][ilat][0] || z < zs[ilon][ilat][nz - 1] ||
      z > zs[ilon][ilat + 1][0] || z < zs[ilon][ilat + 1][nz - 1] ||
      z > zs[ilon + 1][ilat][0] || z < zs[ilon + 1][ilat][nz - 1] ||
      z > zs[ilon + 1][ilat + 1][0] || z < zs[ilon + 1][ilat + 1][nz - 1])
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

  static double xc[NLON][NLAT][3], scal;

  static float helpp[NLON][NLAT], helpt[NLON][NLAT], helpz[NLON][NLAT],
    w, wsum;

  int iz, ilon, ilon2, ilon3, ilat, ilat2, dlon = 3, dlat = 3;

  /* Get Cartesian coordinates... */
  for (ilon = 0; ilon < nlon; ilon++)
    for (ilat = 0; ilat < nlat; ilat++)
      geo2cart(0, lons[ilon], lats[ilat], xc[ilon][ilat]);

  /* Set scaling factor... */
  scal = 1. / (2. * POW2(20. / 2.35482));

  /* Loop over height levels... */
  for (iz = 0; iz < nz; iz++) {

    /* Write info... */
    printf("Smoothing level %d / %d ...\n", iz + 1, nz);

    /* Copy data... */
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++) {
	helpp[ilon][ilat] = ps[ilon][ilat][iz];
	helpt[ilon][ilat] = ts[ilon][ilat][iz];
	helpz[ilon][ilat] = zs[ilon][ilat][iz];
      }

    /* Smoothing... */
    for (ilon = 0; ilon < nlon; ilon++)
      for (ilat = 0; ilat < nlat; ilat++) {
	wsum = 0;
	ps[ilon][ilat][iz] = 0;
	ts[ilon][ilat][iz] = 0;
	zs[ilon][ilat][iz] = 0;
	for (ilon2 = ilon - dlon; ilon2 <= ilon + dlon; ilon2++)
	  for (ilat2 = GSL_MAX(ilat - dlat, 0);
	       ilat2 <= GSL_MIN(ilat + dlat, nlat - 1); ilat2++) {
	    ilon3 = ilon2;
	    if (ilon3 < 0)
	      ilon3 += nlon;
	    else if (ilon3 >= nlon)
	      ilon3 -= nlon;
	    w = (float) exp(-scal * DIST2(xc[ilon][ilat], xc[ilon3][ilat2]));
	    ps[ilon][ilat][iz] += w * helpp[ilon3][ilat2];
	    ts[ilon][ilat][iz] += w * helpt[ilon3][ilat2];
	    zs[ilon][ilat][iz] += w * helpz[ilon3][ilat2];
	    wsum += w;
	  }
	ps[ilon][ilat][iz] /= wsum;
	ts[ilon][ilat][iz] /= wsum;
	zs[ilon][ilat][iz] /= wsum;
      }
  }
}
