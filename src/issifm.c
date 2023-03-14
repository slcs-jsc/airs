#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of model levels. */
#define NZ 152

/*! Maximum number of model longitudes. */
#define NLON 3502

/*! Maximum number of model latitudes. */
#define NLAT 1402

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! Model data. */
typedef struct {

  /*! Number of vertical levels. */
  int nz;

  /*! Number of longitudes. */
  int nlon;

  /*! Number of latitudes. */
  int nlat;

  /*! Longitude [deg]. */
  double lon[NLON];

  /*! Latitude [deg]. */
  double lat[NLAT];

  /*! Surface pressure [hPa]. */
  float ps[NLON][NLAT];

  /*! Pressure [hPa]. */
  float p[NLON][NLAT][NZ];

  /*! Temperature [K]. */
  float t[NLON][NLAT][NZ];

  /*! Height [km]. */
  float z[NLON][NLAT][NZ];

} model_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Interpolation of model data. */
void intpol(
  model_t * model,
  double z,
  double lon,
  double lat,
  double *p,
  double *t);

/*! Smoothing of model data. */
void smooth(
  model_t * model);

/*! Write wave struct to netCDF file. */
void write_nc(
  char *filename,
  wave_t * wave);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;

  static char kernel[LEN], pertname[LEN];

  static double xo[3], xs[3], xm[3], var_dh = 100.,
    f, t_ovp, hyam[NZ], hybm[NZ], kz[NSHAPE], kw[NSHAPE], w, wsum;

  static float *help;

  static int init, id, itrack, ixtrack, ncid, dimid, varid, slant,
    ilon, ilat, iz, ip, track0, track1, nk, okay, extpol;

  static size_t rs;

  atm_t *atm;

  obs_t *obs;

  pert_t *pert;

  wave_t *wave;

  model_t *model;

  /* ------------------------------------------------------------
     Get control parameters...
     ------------------------------------------------------------ */

  /* Check arguments... */
  if (argc < 8)
    ERRMSG("Give parameters: <ctl> <model> <model.nc> <pert.nc>"
	   " <wave_airs.tab> <wave_model.tab> <wave_airs.nc> <wave_model.nc>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  scan_ctl(argc, argv, "KERNEL", -1, "-", kernel);
  extpol = (int) scan_ctl(argc, argv, "EXTPOL", -1, "0", NULL);
  slant = (int) scan_ctl(argc, argv, "SLANT", -1, "1", NULL);
  t_ovp = scan_ctl(argc, argv, "T_OVP", -1, "", NULL);

  /* Set control parameters... */
  ctl.write_bbt = 1;

  /* ------------------------------------------------------------
     Read model data...
     ------------------------------------------------------------ */

  /* Allocate... */
  ALLOC(model, model_t, 1);
  ALLOC(help, float,
	NLON * NLAT * NZ);

  /* Open file... */
  printf("Read %s data: %s\n", argv[2], argv[3]);
  NC(nc_open(argv[3], NC_NOWRITE, &ncid));

  /* Check time... */
  LOG(2, "Check time...");
  if (nc_inq_dimid(ncid, "time", &dimid) != NC_NOERR)
    NC(nc_inq_dimid(ncid, "time", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  if (rs != 1)
    ERRMSG("Only one time step is allowed!");

  /* Read latitudes... */
  LOG(2, "Read latitudes...");
  if (nc_inq_dimid(ncid, "lat", &dimid) != NC_NOERR)
    NC(nc_inq_dimid(ncid, "latitude", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  model->nlat = (int) rs;
  if (model->nlat > NLAT)
    ERRMSG("Too many latitudes!");
  if (nc_inq_varid(ncid, "lat", &varid) != NC_NOERR)
    NC(nc_inq_varid(ncid, "latitude", &varid));
  NC(nc_get_var_double(ncid, varid, model->lat));

  /* Read longitudes... */
  LOG(2, "Read longitudes...");
  if (nc_inq_dimid(ncid, "lon", &dimid) != NC_NOERR)
    NC(nc_inq_dimid(ncid, "longitude", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &rs));
  model->nlon = (int) rs;
  if (model->nlon > NLON)
    ERRMSG("Too many longitudes!");
  if (nc_inq_varid(ncid, "lon", &varid))
    NC(nc_inq_varid(ncid, "longitude", &varid));
  NC(nc_get_var_double(ncid, varid, model->lon));

  /* Read ICON data... */
  if (strcasecmp(argv[2], "icon") == 0) {

    /* Get height levels... */
    LOG(2, "Read levels...");
    NC(nc_inq_dimid(ncid, "height", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    model->nz = (int) rs;
    if (model->nz > NZ)
      ERRMSG("Too many altitudes!");

    /* Read height... */
    LOG(2, "Read height...");
    NC(nc_inq_varid(ncid, "z_mc", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->z[ilon][ilat][iz] =
	    (float) (help[(iz * model->nlat + ilat) * model->nlon + ilon] /
		     1e3);

    /* Read temperature... */
    LOG(2, "Read temperature...");
    NC(nc_inq_varid(ncid, "temp", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->t[ilon][ilat][iz] =
	    help[(iz * model->nlat + ilat) * model->nlon + ilon];

    /* Read pressure... */
    LOG(2, "Read pressure...");
    NC(nc_inq_varid(ncid, "pres", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->p[ilon][ilat][iz] =
	    (float) (help[(iz * model->nlat + ilat) * model->nlon + ilon] /
		     1e2);
  }

  /* Read IFS data... */
  else if (strcasecmp(argv[2], "ifs") == 0) {

    /* Get height levels... */
    LOG(2, "Read levels...");
    NC(nc_inq_dimid(ncid, "lev_2", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    model->nz = (int) rs;
    if (model->nz > NZ)
      ERRMSG("Too many altitudes!");

    /* Read height... */
    LOG(2, "Read height...");
    NC(nc_inq_varid(ncid, "gh", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->z[ilon][ilat][iz] =
	    (float) (help[(iz * model->nlat + ilat) * model->nlon + ilon] /
		     1e3);

    /* Read temperature... */
    LOG(2, "Read temperature...");
    NC(nc_inq_varid(ncid, "t", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->t[ilon][ilat][iz] =
	    help[(iz * model->nlat + ilat) * model->nlon + ilon];

    /* Read surface pressure... */
    LOG(2, "Read surface pressure...");
    NC(nc_inq_varid(ncid, "lnsp", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	model->ps[ilon][ilat] = (float) exp(help[ilat * model->nlon + ilon]);

    /* Read grid coefficients... */
    LOG(2, "Read grid coefficients...");
    NC(nc_inq_varid(ncid, "hyam", &varid));
    NC(nc_get_var_double(ncid, varid, hyam));
    NC(nc_inq_varid(ncid, "hybm", &varid));
    NC(nc_get_var_double(ncid, varid, hybm));

    /* Calculate pressure... */
    LOG(2, "Calculate pressure...");
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->p[ilon][ilat][iz]
	    = (float) ((hyam[iz] + hybm[iz] * model->ps[ilon][ilat]) / 100.);
  }

  /* Read UM data... */
  else if (strcasecmp(argv[2], "um") == 0) {

    /* Get height levels... */
    LOG(2, "Read levels...");
    if (nc_inq_dimid(ncid, "RHO_TOP_eta_rho", &dimid) != NC_NOERR)
      NC(nc_inq_dimid(ncid, "RHO_eta_rho", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    model->nz = (int) rs;
    if (model->nz > NZ)
      ERRMSG("Too many altitudes!");

    /* Read height... */
    LOG(2, "Read height...");
    if (nc_inq_varid(ncid, "STASH_m01s15i102_2", &varid) != NC_NOERR)
      NC(nc_inq_varid(ncid, "STASH_m01s15i102", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->z[ilon][ilat][iz] =
	    (float) (help[(iz * model->nlat + ilat) * model->nlon + ilon] /
		     1e3);

    /* Read temperature... */
    LOG(2, "Read temperature...");
    NC(nc_inq_varid(ncid, "STASH_m01s30i004", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->t[ilon][ilat][iz] =
	    help[(iz * model->nlat + ilat) * model->nlon + ilon];

    /* Read pressure... */
    LOG(2, "Read pressure...");
    NC(nc_inq_varid(ncid, "STASH_m01s00i407", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->p[ilon][ilat][iz] =
	    0.01f * help[(iz * model->nlat + ilat) * model->nlon + ilon];
  }

  /* Read WRF data... */
  else if (strcasecmp(argv[2], "wrf") == 0) {

    /* Get height levels... */
    LOG(2, "Read levels...");
    NC(nc_inq_dimid(ncid, "bottom_top", &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &rs));
    model->nz = (int) rs;
    if (model->nz > NZ)
      ERRMSG("Too many altitudes!");

    /* Read height... */
    LOG(2, "Read height...");
    NC(nc_inq_varid(ncid, "z", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->z[ilon][ilat][iz] =
	    (float) (help[(iz * model->nlat + ilat) * model->nlon + ilon] /
		     1e3);

    /* Read temperature... */
    LOG(2, "Read temperature...");
    NC(nc_inq_varid(ncid, "tk", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->t[ilon][ilat][iz] =
	    help[(iz * model->nlat + ilat) * model->nlon + ilon];

    /* Read pressure... */
    LOG(2, "Read pressure...");
    NC(nc_inq_varid(ncid, "p", &varid));
    NC(nc_get_var_float(ncid, varid, help));
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++)
	for (iz = 0; iz < model->nz; iz++)
	  model->p[ilon][ilat][iz] =
	    (float) (help[(iz * model->nlat + ilat) * model->nlon + ilon] /
		     1e2);
  }

  else
    ERRMSG("Model type not supported!");

  /* Close file... */
  NC(nc_close(ncid));

  /* Free... */
  free(help);

  /* Check data... */
  LOG(2, "Checking...");
  for (ilon = 0; ilon < model->nlon; ilon++)
    for (ilat = 0; ilat < model->nlat; ilat++)
      for (iz = 0; iz < model->nz; iz++)
	if (model->t[ilon][ilat][iz] <= 100
	    || model->t[ilon][ilat][iz] >= 400) {
	  model->p[ilon][ilat][iz] = GSL_NAN;
	  model->t[ilon][ilat][iz] = GSL_NAN;
	  model->z[ilon][ilat][iz] = GSL_NAN;
	}

  /* Smoothing of model data... */
  LOG(2, "Smoothing...");
  smooth(model);

  /* Write info... */
  for (iz = 0; iz < model->nz; iz++)
    printf("section_height: %d %g %g %g %g %g\n", iz,
	   model->z[model->nlon / 2][model->nlat / 2][iz],
	   model->lon[model->nlon / 2], model->lat[model->nlat / 2],
	   model->p[model->nlon / 2][model->nlat / 2][iz],
	   model->t[model->nlon / 2][model->nlat / 2][iz]);
  for (ilon = 0; ilon < model->nlon; ilon++)
    printf("section_west_east: %d %g %g %g %g %g\n", ilon,
	   model->z[ilon][model->nlat / 2][model->nz / 2],
	   model->lon[ilon], model->lat[model->nlat / 2],
	   model->p[ilon][model->nlat / 2][model->nz / 2],
	   model->t[ilon][model->nlat / 2][model->nz / 2]);
  for (ilat = 0; ilat < model->nlat; ilat++)
    printf("section_north_south: %d %g %g %g %g %g\n", ilat,
	   model->z[model->nlon / 2][ilat][model->nz / 2],
	   model->lon[model->nlon / 2], model->lat[ilat],
	   model->p[model->nlon / 2][ilat][model->nz / 2],
	   model->t[model->nlon / 2][ilat][model->nz / 2]);

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
  write_nc(argv[7], wave);

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
      obs->vpz[0] = 0;
      obs->vplon[0] = pert->lon[itrack][ixtrack];
      obs->vplat[0] = pert->lat[itrack][ixtrack];

      /* Get Cartesian coordinates... */
      geo2cart(obs->obsz[0], obs->obslon[0], obs->obslat[0], xo);
      geo2cart(obs->vpz[0], obs->vplon[0], obs->vplat[0], xs);

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
	intpol(model, atm->z[ip], atm->lon[ip], atm->lat[ip],
	       &atm->p[ip], &atm->t[ip]);

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

  /* Extrapolate... */
  if (extpol)
    for (itrack = track0; itrack <= track1; itrack++) {
      for (ixtrack = 1; ixtrack < pert->nxtrack; ixtrack++)
	if (!gsl_finite(pert->bt[itrack][ixtrack]))
	  pert->bt[itrack][ixtrack] = pert->bt[itrack][ixtrack - 1];
      for (ixtrack = pert->nxtrack - 2; ixtrack >= 0; ixtrack--)
	if (!gsl_finite(pert->bt[itrack][ixtrack]))
	  pert->bt[itrack][ixtrack] = pert->bt[itrack][ixtrack + 1];
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
  write_nc(argv[8], wave);

  /* Free... */
  free(atm);
  free(obs);
  free(pert);
  free(wave);

  return EXIT_SUCCESS;
}

/************************************************************************/

void intpol(
  model_t * model,
  double z,
  double lon,
  double lat,
  double *p,
  double *t) {

  double p00, p01, p10, p11, t00, t01, t10, t11, zd[NZ];

  int iz, ilon, ilat;

  /* Adjust longitude... */
  if (model->lon[model->nlon - 1] > 180)
    if (lon < 0)
      lon += 360;

  /* Check horizontal range... */
  if (lon < model->lon[0]
      || lon > model->lon[model->nlon - 1]
      || lat < GSL_MIN(model->lat[0], model->lat[model->nlat - 1])
      || lat > GSL_MAX(model->lat[0], model->lat[model->nlat - 1])) {
    *p = GSL_NAN;
    *t = GSL_NAN;
    return;
  }

  /* Get indices... */
  ilon = locate_irr(model->lon, model->nlon, lon);
  ilat = locate_irr(model->lat, model->nlat, lat);

  /* Check data... */
  if (!gsl_finite(model->z[ilon][ilat][0])
      || !gsl_finite(model->z[ilon][ilat][model->nz - 1])
      || !gsl_finite(model->z[ilon][ilat + 1][model->nz - 1])
      || !gsl_finite(model->z[ilon][ilat + 1][model->nz - 1])
      || !gsl_finite(model->z[ilon + 1][ilat][model->nz - 1])
      || !gsl_finite(model->z[ilon + 1][ilat][model->nz - 1])
      || !gsl_finite(model->z[ilon + 1][ilat + 1][model->nz - 1])
      || !gsl_finite(model->z[ilon + 1][ilat + 1][model->nz - 1])) {
    *p = GSL_NAN;
    *t = GSL_NAN;
    return;
  }

  /* Check vertical range... */
  if (z >
      GSL_MAX(model->z[ilon][ilat][0], model->z[ilon][ilat][model->nz - 1])
      || z < GSL_MIN(model->z[ilon][ilat][0],
		     model->z[ilon][ilat][model->nz - 1])
      || z > GSL_MAX(model->z[ilon][ilat + 1][0],
		     model->z[ilon][ilat + 1][model->nz - 1])
      || z < GSL_MIN(model->z[ilon][ilat + 1][0],
		     model->z[ilon][ilat + 1][model->nz - 1])
      || z > GSL_MAX(model->z[ilon + 1][ilat][0],
		     model->z[ilon + 1][ilat][model->nz - 1])
      || z < GSL_MIN(model->z[ilon + 1][ilat][0],
		     model->z[ilon + 1][ilat][model->nz - 1])
      || z > GSL_MAX(model->z[ilon + 1][ilat + 1][0],
		     model->z[ilon + 1][ilat + 1][model->nz - 1])
      || z < GSL_MIN(model->z[ilon + 1][ilat + 1][0],
		     model->z[ilon + 1][ilat + 1][model->nz - 1]))
    return;

  /* Interpolate vertically... */
  for (iz = 0; iz < model->nz; iz++)
    zd[iz] = model->z[ilon][ilat][iz];
  iz = locate_irr(zd, model->nz, z);
  p00 = LIN(model->z[ilon][ilat][iz], model->p[ilon][ilat][iz],
	    model->z[ilon][ilat][iz + 1], model->p[ilon][ilat][iz + 1], z);
  t00 = LIN(model->z[ilon][ilat][iz], model->t[ilon][ilat][iz],
	    model->z[ilon][ilat][iz + 1], model->t[ilon][ilat][iz + 1], z);

  for (iz = 0; iz < model->nz; iz++)
    zd[iz] = model->z[ilon][ilat + 1][iz];
  iz = locate_irr(zd, model->nz, z);
  p01 = LIN(model->z[ilon][ilat + 1][iz], model->p[ilon][ilat + 1][iz],
	    model->z[ilon][ilat + 1][iz + 1],
	    model->p[ilon][ilat + 1][iz + 1], z);
  t01 =
    LIN(model->z[ilon][ilat + 1][iz], model->t[ilon][ilat + 1][iz],
	model->z[ilon][ilat + 1][iz + 1], model->t[ilon][ilat + 1][iz + 1],
	z);

  for (iz = 0; iz < model->nz; iz++)
    zd[iz] = model->z[ilon + 1][ilat][iz];
  iz = locate_irr(zd, model->nz, z);
  p10 = LIN(model->z[ilon + 1][ilat][iz], model->p[ilon + 1][ilat][iz],
	    model->z[ilon + 1][ilat][iz + 1],
	    model->p[ilon + 1][ilat][iz + 1], z);
  t10 =
    LIN(model->z[ilon + 1][ilat][iz], model->t[ilon + 1][ilat][iz],
	model->z[ilon + 1][ilat][iz + 1], model->t[ilon + 1][ilat][iz + 1],
	z);

  for (iz = 0; iz < model->nz; iz++)
    zd[iz] = model->z[ilon + 1][ilat + 1][iz];
  iz = locate_irr(zd, model->nz, z);
  p11 =
    LIN(model->z[ilon + 1][ilat + 1][iz], model->p[ilon + 1][ilat + 1][iz],
	model->z[ilon + 1][ilat + 1][iz + 1],
	model->p[ilon + 1][ilat + 1][iz + 1], z);
  t11 =
    LIN(model->z[ilon + 1][ilat + 1][iz], model->t[ilon + 1][ilat + 1][iz],
	model->z[ilon + 1][ilat + 1][iz + 1],
	model->t[ilon + 1][ilat + 1][iz + 1], z);

  /* Interpolate horizontally... */
  p00 = LIN(model->lon[ilon], p00, model->lon[ilon + 1], p10, lon);
  p11 = LIN(model->lon[ilon], p01, model->lon[ilon + 1], p11, lon);
  *p = LIN(model->lat[ilat], p00, model->lat[ilat + 1], p11, lat);

  t00 = LIN(model->lon[ilon], t00, model->lon[ilon + 1], t10, lon);
  t11 = LIN(model->lon[ilon], t01, model->lon[ilon + 1], t11, lon);
  *t = LIN(model->lat[ilat], t00, model->lat[ilat + 1], t11, lat);
}

/************************************************************************/

void smooth(
  model_t * model) {

  static float hp[NLON][NLAT], ht[NLON][NLAT], hz[NLON][NLAT], w, wsum;

  static double dx, dy, wx[10], wy[10];

  int iz, ilon, ilon2, ilat, ilat2, dlon = 3, dlat = 3;

  /* Set weights... */
  dy = RE * M_PI / 180. * fabs(model->lat[1] - model->lat[0]);
  for (ilat = 0; ilat <= dlat; ilat++)
    wy[ilat] = exp(-0.5 * POW2(ilat * dy * 2.35482 / 20.));

  /* Loop over height levels... */
  for (iz = 0; iz < model->nz; iz++) {

    /* Write info... */
    printf("Smoothing level %d / %d ...\n", iz + 1, model->nz);

    /* Copy data... */
    for (ilon = 0; ilon < model->nlon; ilon++)
      for (ilat = 0; ilat < model->nlat; ilat++) {
	hp[ilon][ilat] = model->p[ilon][ilat][iz];
	ht[ilon][ilat] = model->t[ilon][ilat][iz];
	hz[ilon][ilat] = model->z[ilon][ilat][iz];
      }

    /* Loop over latitudes... */
    for (ilat = 0; ilat < model->nlat; ilat++) {

      /* Set weights... */
      dx = RE * M_PI / 180. * cos(model->lat[ilat] * M_PI / 180.) *
	fabs(model->lon[1] - model->lon[0]);
      for (ilon = 0; ilon <= dlon; ilon++)
	wx[ilon] = exp(-0.5 * POW2(ilon * dx * 2.35482 / 20.));

      /* Loop over longitudes... */
      for (ilon = 0; ilon < model->nlon; ilon++) {
	wsum = 0;
	model->p[ilon][ilat][iz] = 0;
	model->t[ilon][ilat][iz] = 0;
	model->z[ilon][ilat][iz] = 0;
	for (ilon2 = GSL_MAX(ilon - dlon, 0);
	     ilon2 <= GSL_MIN(ilon + dlon, model->nlon - 1); ilon2++)
	  for (ilat2 = GSL_MAX(ilat - dlat, 0);
	       ilat2 <= GSL_MIN(ilat + dlat, model->nlat - 1); ilat2++) {
	    w = (float) (wx[abs(ilon2 - ilon)] * wy[abs(ilat2 - ilat)]);
	    model->p[ilon][ilat][iz] += w * hp[ilon2][ilat2];
	    model->t[ilon][ilat][iz] += w * ht[ilon2][ilat2];
	    model->z[ilon][ilat][iz] += w * hz[ilon2][ilat2];
	    wsum += w;
	  }
	model->p[ilon][ilat][iz] /= wsum;
	model->t[ilon][ilat][iz] /= wsum;
	model->z[ilon][ilat][iz] /= wsum;
      }
    }
  }
}

/************************************************************************/

void write_nc(
  char *filename,
  wave_t * wave) {

  static double help[WX * WY];

  int ix, iy, ncid, dimid[10], lon_id, lat_id, bt_id, pt_id, var_id;

  /* Create netCDF file... */
  NC(nc_create(filename, NC_CLOBBER, &ncid));

  /* Set dimensions... */
  NC(nc_def_dim(ncid, "NTRACK", (size_t) wave->ny, &dimid[0]));
  NC(nc_def_dim(ncid, "NXTRACK", (size_t) wave->nx, &dimid[1]));

  /* Add variables... */
  NC(nc_def_var(ncid, "lon", NC_DOUBLE, 2, dimid, &lon_id));
  add_att(ncid, lon_id, "deg", "footprint longitude");
  NC(nc_def_var(ncid, "lat", NC_DOUBLE, 2, dimid, &lat_id));
  add_att(ncid, lat_id, "deg", "footprint latitude");
  NC(nc_def_var(ncid, "bt", NC_FLOAT, 2, dimid, &bt_id));
  add_att(ncid, bt_id, "K", "brightness temperature");
  NC(nc_def_var(ncid, "bt_pt", NC_FLOAT, 2, dimid, &pt_id));
  add_att(ncid, pt_id, "K", "brightness temperature perturbation");
  NC(nc_def_var(ncid, "bt_var", NC_FLOAT, 2, dimid, &var_id));
  add_att(ncid, var_id, "K^2", "brightness temperature variance");

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Write data... */
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++)
      help[iy * wave->nx + ix] = wave->lon[ix][iy];
  NC(nc_put_var_double(ncid, lon_id, help));
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++)
      help[iy * wave->nx + ix] = wave->lat[ix][iy];
  NC(nc_put_var_double(ncid, lat_id, help));
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++)
      help[iy * wave->nx + ix] = wave->temp[ix][iy];
  NC(nc_put_var_double(ncid, bt_id, help));
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++)
      help[iy * wave->nx + ix] = wave->pt[ix][iy];
  NC(nc_put_var_double(ncid, pt_id, help));
  for (ix = 0; ix < wave->nx; ix++)
    for (iy = 0; iy < wave->ny; iy++)
      help[iy * wave->nx + ix] = wave->var[ix][iy];
  NC(nc_put_var_double(ncid, var_id, help));

  /* Close file... */
  NC(nc_close(ncid));
}
