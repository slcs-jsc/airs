#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum model dimensions. */
#define NLON 1441
#define NLAT 721
#define NZ 138

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

void intpol(
  float xs[NLON][NLAT][NZ],
  double zs[NLON][NLAT][NZ],
  double lons[NLON],
  double lats[NLAT],
  int nz,
  int nlon,
  int nlat,
  double z,
  double lon,
  double lat,
  double *x);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static atm_t atm;
  static obs_t obs;

  static char pertname[LEN];

  static double lon[NLON], lat[NLAT], xo[3], xs[3], xm[3], var_dh = 100.,
    f, z[NLON][NLAT][NZ], t_ovp, hyam[NZ], hybm[NZ], ps[NLON][NLAT];

  static float p[NLON][NLAT][NZ], t[NLON][NLAT][NZ], help[NLON * NLAT * NZ];

  static int id, itrack, ixtrack, ncid, dimid, varid,
    ilon, ilat, iz, nlon, nlat, nz, ip, track0, track1;

  static size_t rs;

  pert_t *pert;

  wave_t *wave;

  /* ------------------------------------------------------------
     Get control parameters...
     ------------------------------------------------------------ */

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <model.nc> <pert.nc>"
	   " <wave_airs.tab> <wave_model.tab>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  t_ovp = scan_ctl(argc, argv, "T_OVP", -1, "", NULL);

  /* Set control parameters... */
  ctl.write_bbt = 1;

  /* Allocate... */
  ALLOC(pert, pert_t, 1);
  ALLOC(wave, wave_t, 1);

  /* ------------------------------------------------------------
     Read model data...
     ------------------------------------------------------------ */

  /* Open file... */
  printf("Read model data: %s\n", argv[2]);
  NC(nc_open(argv[2], NC_NOWRITE, &ncid));

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
	z[ilon][ilat][iz] = help[(iz * nlat + ilat) * nlon + ilon] / 1e3;

  /* Read surface pressure... */
  NC(nc_inq_varid(ncid, "lnsp", &varid));
  NC(nc_get_var_float(ncid, varid, help));
  for (ilon = 0; ilon < nlon; ilon++)
    for (ilat = 0; ilat < nlat; ilat++)
      ps[ilon][ilat] = exp(help[ilat * nlon + ilon]);

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

  /* ------------------------------------------------------------
     Read AIRS perturbation data...
     ------------------------------------------------------------ */

  /* Read perturbation data... */
  read_pert(argv[3], pertname, pert);

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
  write_wave(argv[4], wave);

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
      obs.nr = 1;
      obs.obsz[0] = 705;
      obs.obslon[0] = pert->lon[itrack][44];
      obs.obslat[0] = pert->lat[itrack][44];

      /* Get Cartesian coordinates... */
      geo2cart(obs.obsz[0], obs.obslon[0], obs.obslat[0], xo);
      geo2cart(0, pert->lon[itrack][ixtrack], pert->lat[itrack][ixtrack], xs);

      /* Set atmospheric data... */
      atm.np = 0;
      for (f = 0.0; f <= 1.0; f += 0.0002) {
	xm[0] = f * xo[0] + (1 - f) * xs[0];
	xm[1] = f * xo[1] + (1 - f) * xs[1];
	xm[2] = f * xo[2] + (1 - f) * xs[2];
	cart2geo(xm, &atm.z[atm.np], &atm.lon[atm.np], &atm.lat[atm.np]);
	atm.time[atm.np] = pert->time[itrack][ixtrack];
	if (atm.z[atm.np] < 10)
	  continue;
	else if (atm.z[atm.np] > 90)
	  break;
	else if ((++atm.np) >= NP)
	  ERRMSG("Too many altitudes!");
      }

      /* Initialize with climatological data... */
      climatology(&ctl, &atm);

      /* Interpolate model data... */
      for (ip = 0; ip < atm.np; ip++) {
	intpol(t, z, lon, lat, nz, nlon, nlat,
	       atm.z[ip], atm.lon[ip], atm.lat[ip], &atm.t[ip]);
	intpol(p, z, lon, lat, nz, nlon, nlat,
	       atm.z[ip], atm.lon[ip], atm.lat[ip], &atm.p[ip]);
      }

      /* Run forward model... */
      formod(&ctl, &atm, &obs);

      /* Get mean brightness temperature... */
      pert->bt[itrack][ixtrack] = 0;
      for (id = 0; id < ctl.nd; id++)
	pert->bt[itrack][ixtrack] += obs.rad[id][0] / ctl.nd;
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
  write_wave(argv[5], wave);

  /* Free... */
  free(pert);
  free(wave);

  return EXIT_SUCCESS;
}

/************************************************************************/

void intpol(
  float xs[NLON][NLAT][NZ],
  double zs[NLON][NLAT][NZ],
  double lons[NLON],
  double lats[NLAT],
  int nz,
  int nlon,
  int nlat,
  double z,
  double lon,
  double lat,
  double *x) {

  double x00, x01, x10, x11;

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
  iz = locate_irr(zs[ilon][ilat], nz, z);
  x00 = LIN(zs[ilon][ilat][iz], xs[ilon][ilat][iz],
	    zs[ilon][ilat][iz + 1], xs[ilon][ilat][iz + 1], z);

  iz = locate_irr(zs[ilon][ilat + 1], nz, z);
  x01 = LIN(zs[ilon][ilat + 1][iz], xs[ilon][ilat + 1][iz],
	    zs[ilon][ilat + 1][iz + 1], xs[ilon][ilat + 1][iz + 1], z);

  iz = locate_irr(zs[ilon + 1][ilat], nz, z);
  x10 = LIN(zs[ilon + 1][ilat][iz], xs[ilon + 1][ilat][iz],
	    zs[ilon + 1][ilat][iz + 1], xs[ilon + 1][ilat][iz + 1], z);

  iz = locate_irr(zs[ilon + 1][ilat + 1], nz, z);
  x11 = LIN(zs[ilon + 1][ilat + 1][iz], xs[ilon + 1][ilat + 1][iz],
	    zs[ilon + 1][ilat + 1][iz + 1], xs[ilon + 1][ilat + 1][iz + 1],
	    z);

  /* Interpolate horizontally... */
  x00 = LIN(lons[ilon], x00, lons[ilon + 1], x10, lon);
  x11 = LIN(lons[ilon], x01, lons[ilon + 1], x11, lon);
  *x = LIN(lats[ilat], x00, lats[ilat + 1], x11, lat);
}
