#include "libairs.h"

/* ------------------------------------------------------------
   Global variables...
   ------------------------------------------------------------ */

/* List of AIRS channels (don't change). */
int airs_chan[L1_NCHAN] = { 54, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
  2035, 2036, 2040, 2041, 2052, 2053, 2054, 2055,
  2067, 2075, 2076, 2077, 2078, 2079, 2080, 2081,
  2082, 2086, 2088, 2089, 2091, 2092, 2093
};

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of pressure levels for meteorological data. */
#define EP 72

/* Maximum number of longitudes for meteorological data. */
#define EX 362

/* Maximum number of latitudes for meteorological data. */
#define EY 182

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/* Control parameters. */
typedef struct {

  /* Time step of meteorological data [s]. */
  double dt_met;

  /* Surface geopotential data file. */
  char met_geopot[LEN];

} ctl2_t;

/* Meteorological data. */
typedef struct {

  /* Time [s]. */
  double time;

  /* Number of longitudes. */
  int nx;

  /* Number of latitudes. */
  int ny;

  /* Number of pressure levels. */
  int np;

  /* Longitude [deg]. */
  double lon[EX];

  /* Latitude [deg]. */
  double lat[EY];

  /* Pressure [hPa]. */
  double p[EP];

  /* Surface pressure [hPa]. */
  double ps[EX][EY];

  /* Tropopause pressure [hPa]. */
  double pt[EX][EY];

  /* Geopotential height [km]. */
  float z[EX][EY][EP];

  /* Temperature [K]. */
  float t[EX][EY][EP];

  /* Zonal wind [m/s]. */
  float u[EX][EY][EP];

  /* Meridional wind [m/s]. */
  float v[EX][EY][EP];

  /* Vertical wind [hPa/s]. */
  float w[EX][EY][EP];

  /* Potential vorticity [PVU]. */
  float pv[EX][EY][EP];

  /* Water vapor volume mixing ratio [1]. */
  float h2o[EX][EY][EP];

  /* Ozone volume mixing ratio [1]. */
  float o3[EX][EY][EP];

  /* Pressure on model levels [hPa]. */
  float pl[EX][EY][EP];

} met_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Get meteorological data for given timestep. */
void get_met(
  ctl2_t * ctl2,
  char *metbase,
  double t,
  met_t * met0,
  met_t * met1);

/* Get meteorological data for timestep. */
void get_met_help(
  double t,
  int direct,
  char *metbase,
  double dt_met,
  char *filename);

/* Linear interpolation of 2-D meteorological data. */
void intpol_met_2d(
  double array[EX][EY],
  int ix,
  int iy,
  double wx,
  double wy,
  double *var);

/* Linear interpolation of 3-D meteorological data. */
void intpol_met_3d(
  float array[EX][EY][EP],
  int ip,
  int ix,
  int iy,
  double wp,
  double wx,
  double wy,
  double *var);

/* Spatial interpolation of meteorological data. */
void intpol_met_space(
  met_t * met,
  double p,
  double lon,
  double lat,
  double *ps,
  double *pt,
  double *z,
  double *t,
  double *u,
  double *v,
  double *w,
  double *pv,
  double *h2o,
  double *o3);

/* Temporal interpolation of meteorological data. */
void intpol_met_time(
  met_t * met0,
  met_t * met1,
  double ts,
  double p,
  double lon,
  double lat,
  double *ps,
  double *pt,
  double *z,
  double *t,
  double *u,
  double *v,
  double *w,
  double *pv,
  double *h2o,
  double *o3);

/* Read control parameters. */
void read_ctl2(
  int argc,
  char *argv[],
  ctl2_t * ctl2);

/* Read meteorological data file. */
void read_met(
  ctl2_t * ctl2,
  char *filename,
  met_t * met);

/* Extrapolate meteorological data at lower boundary. */
void read_met_extrapolate(
  met_t * met);

/* Calculate geopotential heights. */
void read_met_geopot(
  ctl2_t * ctl2,
  met_t * met);

/* Read and convert variable from meteorological data file. */
void read_met_help(
  int ncid,
  char *varname,
  char *varname2,
  met_t * met,
  float dest[EX][EY][EP],
  float scl);

/* Create meteorological data with periodic boundary conditions. */
void read_met_periodic(
  met_t * met);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;

  static airs_l1_t l1;
  static airs_l2_t l2;

  static ctl2_t ctl2;

  met_t *met0, *met1;

  double ts;

  int ichan, lay, track, xtrack;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <airs_l1_file> <metbase> <out.nc>");

  /* Allocate... */
  ALLOC(met0, met_t, 1);
  ALLOC(met1, met_t, 1);

  /* Read control parameters... */
  read_ctl2(argc, argv, &ctl2);

  /* Read AIRS data... */
  printf("Read AIRS Level-1 file: %s\n", argv[2]);
  airs_rad_rdr(argv[2], &airs_rad_gran);

  /* Flag bad data... */
  for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
    for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++)
      for (ichan = 0; ichan < L1_NCHAN; ichan++)
	if ((airs_rad_gran.state[track][xtrack] != 0)
	    || (airs_rad_gran.ExcludedChans[airs_chan[ichan]] > 2)
	    || (airs_rad_gran.CalChanSummary[airs_chan[ichan]] & 8)
	    || (airs_rad_gran.CalChanSummary[airs_chan[ichan]] & (32 + 64))
	    || (airs_rad_gran.CalFlag[track][airs_chan[ichan]] & 16))
	  airs_rad_gran.radiances[track][xtrack][airs_chan[ichan]]
	    = GSL_NAN;

  /* Copy data to struct... */
  for (track = 0; track < AIRS_RAD_GEOTRACK; track++)
    for (xtrack = 0; xtrack < AIRS_RAD_GEOXTRACK; xtrack++) {
      l1.time[track][xtrack]
	= airs_rad_gran.Time[track][xtrack] - 220838400.;
      l1.lon[track][xtrack]
	= airs_rad_gran.Longitude[track][xtrack];
      l1.lat[track][xtrack]
	= airs_rad_gran.Latitude[track][xtrack];
      l1.sat_z[track]
	= airs_rad_gran.satheight[track];
      l1.sat_lon[track]
	= airs_rad_gran.sat_lon[track];
      l1.sat_lat[track]
	= airs_rad_gran.sat_lat[track];
      for (ichan = 0; ichan < L1_NCHAN; ichan++) {
	l1.nu[ichan]
	  = airs_rad_gran.nominal_freq[airs_chan[ichan]];
	l1.rad[track][xtrack][ichan]
	  = airs_rad_gran.radiances[track][xtrack][airs_chan[ichan]] * 0.001f;
      }
    }

  /* Write netCDF file... */
  write_l1(argv[4], &l1);

  /* Read meteo data... */
  ts =
    airs_rad_gran.Time[AIRS_RAD_GEOTRACK / 2][AIRS_RAD_GEOXTRACK / 2] -
    220838400.;
  get_met(&ctl2, argv[3], ts, met0, met1);

  /* Interpolate meteo data... */
  for (track = 0; track < L2_NTRACK; track++)
    for (xtrack = 0; xtrack < L2_NXTRACK; xtrack++)
      for (lay = 0; lay < L2_NLAY; lay++) {
	l2.time[track][xtrack]
	  = airs_rad_gran.Time[track * 3 + 1][xtrack * 3 + 1] - 220838400.;
	l2.lon[track][xtrack]
	  = airs_rad_gran.Longitude[track * 3 + 1][xtrack * 3 + 1];
	l2.lat[track][xtrack]
	  = airs_rad_gran.Latitude[track * 3 + 1][xtrack * 3 + 1];
	l2.p[lay]
	  = 1013.25 * exp(-2.5 * lay / 7.0);
	intpol_met_time(met0, met1, ts, l2.p[lay],
			l2.lon[track][xtrack], l2.lat[track][xtrack],
			NULL, NULL, &l2.z[track][xtrack][lay],
			&l2.t[track][xtrack][lay],
			NULL, NULL, NULL, NULL, NULL, NULL);
      }

  /* Write netCDF file... */
  write_l2(argv[4], &l2);

  /* Free... */
  free(met0);
  free(met1);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void get_met(
  ctl2_t * ctl2,
  char *metbase,
  double t,
  met_t * met0,
  met_t * met1) {

  char filename[LEN];

  static int init;

  /* Init... */
  if (!init) {
    init = 1;

    get_met_help(t, -1, metbase, ctl2->dt_met, filename);
    read_met(ctl2, filename, met0);

    get_met_help(t + 1.0, 1, metbase, ctl2->dt_met, filename);
    read_met(ctl2, filename, met1);
  }

  /* Read new data for forward trajectories... */
  if (t > met1->time) {
    memcpy(met0, met1, sizeof(met_t));
    get_met_help(t, 1, metbase, ctl2->dt_met, filename);
    read_met(ctl2, filename, met1);
  }
}

/*****************************************************************************/

void get_met_help(
  double t,
  int direct,
  char *metbase,
  double dt_met,
  char *filename) {

  double t6, r;

  int year, mon, day, hour, min, sec;

  /* Round time to fixed intervals... */
  if (direct == -1)
    t6 = floor(t / dt_met) * dt_met;
  else
    t6 = ceil(t / dt_met) * dt_met;

  /* Decode time... */
  jsec2time(t6, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Set filename... */
  sprintf(filename, "%s_%d_%02d_%02d_%02d.nc", metbase, year, mon, day, hour);
}

/*****************************************************************************/

void intpol_met_2d(
  double array[EX][EY],
  int ix,
  int iy,
  double wx,
  double wy,
  double *var) {

  double aux00, aux01, aux10, aux11;

  /* Set variables... */
  aux00 = array[ix][iy];
  aux01 = array[ix][iy + 1];
  aux10 = array[ix + 1][iy];
  aux11 = array[ix + 1][iy + 1];

  /* Interpolate horizontally... */
  aux00 = wy * (aux00 - aux01) + aux01;
  aux11 = wy * (aux10 - aux11) + aux11;
  *var = wx * (aux00 - aux11) + aux11;
}

/*****************************************************************************/

void intpol_met_3d(
  float array[EX][EY][EP],
  int ip,
  int ix,
  int iy,
  double wp,
  double wx,
  double wy,
  double *var) {

  double aux00, aux01, aux10, aux11;

  /* Interpolate vertically... */
  aux00 = wp * (array[ix][iy][ip] - array[ix][iy][ip + 1])
    + array[ix][iy][ip + 1];
  aux01 = wp * (array[ix][iy + 1][ip] - array[ix][iy + 1][ip + 1])
    + array[ix][iy + 1][ip + 1];
  aux10 = wp * (array[ix + 1][iy][ip] - array[ix + 1][iy][ip + 1])
    + array[ix + 1][iy][ip + 1];
  aux11 = wp * (array[ix + 1][iy + 1][ip] - array[ix + 1][iy + 1][ip + 1])
    + array[ix + 1][iy + 1][ip + 1];

  /* Interpolate horizontally... */
  aux00 = wy * (aux00 - aux01) + aux01;
  aux11 = wy * (aux10 - aux11) + aux11;
  *var = wx * (aux00 - aux11) + aux11;
}

/*****************************************************************************/

void intpol_met_space(
  met_t * met,
  double p,
  double lon,
  double lat,
  double *ps,
  double *pt,
  double *z,
  double *t,
  double *u,
  double *v,
  double *w,
  double *pv,
  double *h2o,
  double *o3) {

  double wp, wx, wy;

  int ip, ix, iy;

  /* Check longitude... */
  if (met->lon[met->nx - 1] > 180 && lon < 0)
    lon += 360;

  /* Get indices... */
  ip = locate(met->p, met->np, p);
  ix = locate(met->lon, met->nx, lon);
  iy = locate(met->lat, met->ny, lat);

  /* Get weights... */
  wp = (met->p[ip + 1] - p) / (met->p[ip + 1] - met->p[ip]);
  wx = (met->lon[ix + 1] - lon) / (met->lon[ix + 1] - met->lon[ix]);
  wy = (met->lat[iy + 1] - lat) / (met->lat[iy + 1] - met->lat[iy]);

  /* Interpolate... */
  if (ps != NULL)
    intpol_met_2d(met->ps, ix, iy, wx, wy, ps);
  if (pt != NULL)
    intpol_met_2d(met->pt, ix, iy, wx, wy, pt);
  if (z != NULL)
    intpol_met_3d(met->z, ip, ix, iy, wp, wx, wy, z);
  if (t != NULL)
    intpol_met_3d(met->t, ip, ix, iy, wp, wx, wy, t);
  if (u != NULL)
    intpol_met_3d(met->u, ip, ix, iy, wp, wx, wy, u);
  if (v != NULL)
    intpol_met_3d(met->v, ip, ix, iy, wp, wx, wy, v);
  if (w != NULL)
    intpol_met_3d(met->w, ip, ix, iy, wp, wx, wy, w);
  if (pv != NULL)
    intpol_met_3d(met->pv, ip, ix, iy, wp, wx, wy, pv);
  if (h2o != NULL)
    intpol_met_3d(met->h2o, ip, ix, iy, wp, wx, wy, h2o);
  if (o3 != NULL)
    intpol_met_3d(met->o3, ip, ix, iy, wp, wx, wy, o3);
}

/*****************************************************************************/

void intpol_met_time(
  met_t * met0,
  met_t * met1,
  double ts,
  double p,
  double lon,
  double lat,
  double *ps,
  double *pt,
  double *z,
  double *t,
  double *u,
  double *v,
  double *w,
  double *pv,
  double *h2o,
  double *o3) {

  double h2o0, h2o1, o30, o31, ps0, ps1, pt0, pt1, pv0, pv1, t0, t1, u0, u1,
    v0, v1, w0, w1, wt, z0, z1;

  /* Spatial interpolation... */
  intpol_met_space(met0, p, lon, lat,
		   ps == NULL ? NULL : &ps0,
		   pt == NULL ? NULL : &pt0,
		   z == NULL ? NULL : &z0,
		   t == NULL ? NULL : &t0,
		   u == NULL ? NULL : &u0,
		   v == NULL ? NULL : &v0,
		   w == NULL ? NULL : &w0,
		   pv == NULL ? NULL : &pv0,
		   h2o == NULL ? NULL : &h2o0, o3 == NULL ? NULL : &o30);
  intpol_met_space(met1, p, lon, lat,
		   ps == NULL ? NULL : &ps1,
		   pt == NULL ? NULL : &pt1,
		   z == NULL ? NULL : &z1,
		   t == NULL ? NULL : &t1,
		   u == NULL ? NULL : &u1,
		   v == NULL ? NULL : &v1,
		   w == NULL ? NULL : &w1,
		   pv == NULL ? NULL : &pv1,
		   h2o == NULL ? NULL : &h2o1, o3 == NULL ? NULL : &o31);

  /* Get weighting factor... */
  wt = (met1->time - ts) / (met1->time - met0->time);

  /* Interpolate... */
  if (ps != NULL)
    *ps = wt * (ps0 - ps1) + ps1;
  if (pt != NULL)
    *pt = wt * (pt0 - pt1) + pt1;
  if (z != NULL)
    *z = wt * (z0 - z1) + z1;
  if (t != NULL)
    *t = wt * (t0 - t1) + t1;
  if (u != NULL)
    *u = wt * (u0 - u1) + u1;
  if (v != NULL)
    *v = wt * (v0 - v1) + v1;
  if (w != NULL)
    *w = wt * (w0 - w1) + w1;
  if (pv != NULL)
    *pv = wt * (pv0 - pv1) + pv1;
  if (h2o != NULL)
    *h2o = wt * (h2o0 - h2o1) + h2o1;
  if (o3 != NULL)
    *o3 = wt * (o30 - o31) + o31;
}

/*****************************************************************************/

void read_ctl2(
  int argc,
  char *argv[],
  ctl2_t * ctl2) {

  /* Meteorological data... */
  ctl2->dt_met = scan_ctl(argc, argv, "DT_MET", -1, "21600", NULL);
  scan_ctl(argc, argv, "MET_GEOPOT", -1, "", ctl2->met_geopot);
}

/*****************************************************************************/

void read_met(
  ctl2_t * ctl2,
  char *filename,
  met_t * met) {

  char levname[LEN], tstr[10];

  static float help[EX * EY];

  int ix, iy, ip, dimid, ncid, varid, year, mon, day, hour;

  size_t np, nx, ny;

  /* Write info... */
  printf("Read meteorological data: %s\n", filename);

  /* Get time from filename... */
  sprintf(tstr, "%.4s", &filename[strlen(filename) - 16]);
  year = atoi(tstr);
  sprintf(tstr, "%.2s", &filename[strlen(filename) - 11]);
  mon = atoi(tstr);
  sprintf(tstr, "%.2s", &filename[strlen(filename) - 8]);
  day = atoi(tstr);
  sprintf(tstr, "%.2s", &filename[strlen(filename) - 5]);
  hour = atoi(tstr);
  time2jsec(year, mon, day, hour, 0, 0, 0, &met->time);

  /* Open netCDF file... */
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "lon", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &nx));
  if (nx < 2 || nx > EX)
    ERRMSG("Number of longitudes out of range!");

  NC(nc_inq_dimid(ncid, "lat", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &ny));
  if (ny < 2 || ny > EY)
    ERRMSG("Number of latitudes out of range!");

  sprintf(levname, "lev");
  NC(nc_inq_dimid(ncid, levname, &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &np));
  if (np == 1) {
    sprintf(levname, "lev_2");
    NC(nc_inq_dimid(ncid, levname, &dimid));
    NC(nc_inq_dimlen(ncid, dimid, &np));
  }
  if (np < 2 || np > EP)
    ERRMSG("Number of levels out of range!");

  /* Store dimensions... */
  met->np = (int) np;
  met->nx = (int) nx;
  met->ny = (int) ny;

  /* Get horizontal grid... */
  NC(nc_inq_varid(ncid, "lon", &varid));
  NC(nc_get_var_double(ncid, varid, met->lon));
  NC(nc_inq_varid(ncid, "lat", &varid));
  NC(nc_get_var_double(ncid, varid, met->lat));

  /* Read meteorological data... */
  read_met_help(ncid, "t", "T", met, met->t, 1.0);
  read_met_help(ncid, "u", "U", met, met->u, 1.0);
  read_met_help(ncid, "v", "V", met, met->v, 1.0);
  read_met_help(ncid, "w", "W", met, met->w, 0.01f);
  read_met_help(ncid, "q", "Q", met, met->h2o, 1.608f);
  read_met_help(ncid, "o3", "O3", met, met->o3, 0.602f);

  /* Read pressure levels from file... */
  NC(nc_inq_varid(ncid, levname, &varid));
  NC(nc_get_var_double(ncid, varid, met->p));
  for (ip = 0; ip < met->np; ip++)
    met->p[ip] /= 100.;

  /* Extrapolate data for lower boundary... */
  read_met_extrapolate(met);

  /* Check ordering of pressure levels... */
  for (ip = 1; ip < met->np; ip++)
    if (met->p[ip - 1] < met->p[ip])
      ERRMSG("Pressure levels must be descending!");

  /* Read surface pressure... */
  if (nc_inq_varid(ncid, "ps", &varid) == NC_NOERR
      || nc_inq_varid(ncid, "PS", &varid) == NC_NOERR) {
    NC(nc_get_var_float(ncid, varid, help));
    for (iy = 0; iy < met->ny; iy++)
      for (ix = 0; ix < met->nx; ix++)
	met->ps[ix][iy] = help[iy * met->nx + ix] / 100.;
  } else if (nc_inq_varid(ncid, "lnsp", &varid) == NC_NOERR
	     || nc_inq_varid(ncid, "LNSP", &varid) == NC_NOERR) {
    NC(nc_get_var_float(ncid, varid, help));
    for (iy = 0; iy < met->ny; iy++)
      for (ix = 0; ix < met->nx; ix++)
	met->ps[ix][iy] = exp(help[iy * met->nx + ix]) / 100.;
  } else
    for (ix = 0; ix < met->nx; ix++)
      for (iy = 0; iy < met->ny; iy++)
	met->ps[ix][iy] = met->p[0];

  /* Create periodic boundary conditions... */
  read_met_periodic(met);

  /* Calculate geopotential heights... */
  read_met_geopot(ctl2, met);

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_met_extrapolate(
  met_t * met) {

  int ip, ip0, ix, iy;

  /* Loop over columns... */
  for (ix = 0; ix < met->nx; ix++)
    for (iy = 0; iy < met->ny; iy++) {

      /* Find lowest valid data point... */
      for (ip0 = met->np - 1; ip0 >= 0; ip0--)
	if (!gsl_finite(met->t[ix][iy][ip0])
	    || !gsl_finite(met->u[ix][iy][ip0])
	    || !gsl_finite(met->v[ix][iy][ip0])
	    || !gsl_finite(met->w[ix][iy][ip0]))
	  break;

      /* Extrapolate... */
      for (ip = ip0; ip >= 0; ip--) {
	met->t[ix][iy][ip] = met->t[ix][iy][ip + 1];
	met->u[ix][iy][ip] = met->u[ix][iy][ip + 1];
	met->v[ix][iy][ip] = met->v[ix][iy][ip + 1];
	met->w[ix][iy][ip] = met->w[ix][iy][ip + 1];
	met->h2o[ix][iy][ip] = met->h2o[ix][iy][ip + 1];
	met->o3[ix][iy][ip] = met->o3[ix][iy][ip + 1];
      }
    }
}

/*****************************************************************************/

void read_met_geopot(
  ctl2_t * ctl2,
  met_t * met) {

  static double topo_lat[EY], topo_lon[EX], topo_z[EX][EY];

  static int init, topo_nx = -1, topo_ny;

  FILE *in;

  char line[LEN];

  double data[30], lat, lon, rlat, rlon, rlon_old = -999, rz, ts, z0, z1;

  float help[EX][EY];

  int ip, ip0, ix, ix2, ix3, iy, iy2, n, tx, ty;

  /* Initialize geopotential heights... */
  for (ix = 0; ix < met->nx; ix++)
    for (iy = 0; iy < met->ny; iy++)
      for (ip = 0; ip < met->np; ip++)
	met->z[ix][iy][ip] = GSL_NAN;

  /* Check filename... */
  if (ctl2->met_geopot[0] == '-')
    return;

  /* Read surface geopotential... */
  if (!init) {

    /* Write info... */
    printf("Read surface geopotential: %s\n", ctl2->met_geopot);

    /* Open file... */
    if (!(in = fopen(ctl2->met_geopot, "r")))
      ERRMSG("Cannot open file!");

    /* Read data... */
    while (fgets(line, LEN, in))
      if (sscanf(line, "%lg %lg %lg", &rlon, &rlat, &rz) == 3) {
	if (rlon != rlon_old) {
	  if ((++topo_nx) >= EX)
	    ERRMSG("Too many longitudes!");
	  topo_ny = 0;
	}
	rlon_old = rlon;
	topo_lon[topo_nx] = rlon;
	topo_lat[topo_ny] = rlat;
	topo_z[topo_nx][topo_ny] = rz;
	if ((++topo_ny) >= EY)
	  ERRMSG("Too many latitudes!");
      }
    if ((++topo_nx) >= EX)
      ERRMSG("Too many longitudes!");

    /* Close file... */
    fclose(in);

    /* Check grid spacing... */
    if (fabs(met->lon[0] - met->lon[1]) != fabs(topo_lon[0] - topo_lon[1])
	|| fabs(met->lat[0] - met->lat[1]) != fabs(topo_lat[0] - topo_lat[1]))
      printf("Warning: Grid spacing does not match!\n");

    /* Set init flag... */
    init = 1;
  }

  /* Apply hydrostatic equation to calculate geopotential heights... */
  for (ix = 0; ix < met->nx; ix++)
    for (iy = 0; iy < met->ny; iy++) {

      /* Get surface height... */
      lon = met->lon[ix];
      if (lon < topo_lon[0])
	lon += 360;
      else if (lon > topo_lon[topo_nx - 1])
	lon -= 360;
      lat = met->lat[iy];
      tx = locate(topo_lon, topo_nx, lon);
      ty = locate(topo_lat, topo_ny, lat);
      z0 = LIN(topo_lon[tx], topo_z[tx][ty],
	       topo_lon[tx + 1], topo_z[tx + 1][ty], lon);
      z1 = LIN(topo_lon[tx], topo_z[tx][ty + 1],
	       topo_lon[tx + 1], topo_z[tx + 1][ty + 1], lon);
      z0 = LIN(topo_lat[ty], z0, topo_lat[ty + 1], z1, lat);

      /* Find surface pressure level... */
      ip0 = locate(met->p, met->np, met->ps[ix][iy]);

      /* Get surface temperature... */
      ts = LIN(met->p[ip0], met->t[ix][iy][ip0],
	       met->p[ip0 + 1], met->t[ix][iy][ip0 + 1], met->ps[ix][iy]);

      /* Upper part of profile... */
      met->z[ix][iy][ip0 + 1]
	= (float) (z0 + 8.31441 / 28.9647 / G0
		   * 0.5 * (ts + met->t[ix][iy][ip0 + 1])
		   * log(met->ps[ix][iy] / met->p[ip0 + 1]));
      for (ip = ip0 + 2; ip < met->np; ip++)
	met->z[ix][iy][ip]
	  = (float) (met->z[ix][iy][ip - 1] + 8.31441 / 28.9647 / G0
		     * 0.5 * (met->t[ix][iy][ip - 1] + met->t[ix][iy][ip])
		     * log(met->p[ip - 1] / met->p[ip]));
    }

  /* Smooth fields... */
  for (ip = 0; ip < met->np; ip++) {

    /* Median filter... */
    for (ix = 0; ix < met->nx; ix++)
      for (iy = 0; iy < met->nx; iy++) {
	n = 0;
	for (ix2 = ix - 2; ix2 <= ix + 2; ix2++) {
	  ix3 = ix2;
	  if (ix3 < 0)
	    ix3 += met->nx;
	  if (ix3 >= met->nx)
	    ix3 -= met->nx;
	  for (iy2 = GSL_MAX(iy - 2, 0); iy2 <= GSL_MIN(iy + 2, met->ny - 1);
	       iy2++)
	    if (gsl_finite(met->z[ix3][iy2][ip])) {
	      data[n] = met->z[ix3][iy2][ip];
	      n++;
	    }
	}
	if (n > 0) {
	  gsl_sort(data, 1, (size_t) n);
	  help[ix][iy] = (float)
	    gsl_stats_median_from_sorted_data(data, 1, (size_t) n);
	} else
	  help[ix][iy] = GSL_NAN;
      }

    /* Copy data... */
    for (ix = 0; ix < met->nx; ix++)
      for (iy = 0; iy < met->nx; iy++)
	met->z[ix][iy][ip] = help[ix][iy];
  }
}

/*****************************************************************************/

void read_met_help(
  int ncid,
  char *varname,
  char *varname2,
  met_t * met,
  float dest[EX][EY][EP],
  float scl) {

  static float help[EX * EY * EP];

  int ip, ix, iy, varid;

  /* Check if variable exists... */
  if (nc_inq_varid(ncid, varname, &varid) != NC_NOERR)
    if (nc_inq_varid(ncid, varname2, &varid) != NC_NOERR)
      return;

  /* Read data... */
  NC(nc_get_var_float(ncid, varid, help));

  /* Copy and check data... */
  for (ix = 0; ix < met->nx; ix++)
    for (iy = 0; iy < met->ny; iy++)
      for (ip = 0; ip < met->np; ip++) {
	dest[ix][iy][ip] = help[(ip * met->ny + iy) * met->nx + ix];
	if (fabsf(dest[ix][iy][ip]) < 1e14f)
	  dest[ix][iy][ip] *= scl;
	else
	  dest[ix][iy][ip] = GSL_NAN;
      }
}

/*****************************************************************************/

void read_met_periodic(
  met_t * met) {

  int ip, iy;

  /* Check longitudes... */
  if (!(fabs(met->lon[met->nx - 1] - met->lon[0]
	     + met->lon[1] - met->lon[0] - 360) < 0.01))
    return;

  /* Increase longitude counter... */
  if ((++met->nx) > EX)
    ERRMSG("Cannot create periodic boundary conditions!");

  /* Set longitude... */
  met->lon[met->nx - 1] = met->lon[met->nx - 2] + met->lon[1] - met->lon[0];

  /* Loop over latitudes and pressure levels... */
  for (iy = 0; iy < met->ny; iy++)
    for (ip = 0; ip < met->np; ip++) {
      met->ps[met->nx - 1][iy] = met->ps[0][iy];
      met->pt[met->nx - 1][iy] = met->pt[0][iy];
      met->z[met->nx - 1][iy][ip] = met->z[0][iy][ip];
      met->t[met->nx - 1][iy][ip] = met->t[0][iy][ip];
      met->u[met->nx - 1][iy][ip] = met->u[0][iy][ip];
      met->v[met->nx - 1][iy][ip] = met->v[0][iy][ip];
      met->w[met->nx - 1][iy][ip] = met->w[0][iy][ip];
      met->pv[met->nx - 1][iy][ip] = met->pv[0][iy][ip];
      met->h2o[met->nx - 1][iy][ip] = met->h2o[0][iy][ip];
      met->o3[met->nx - 1][iy][ip] = met->o3[0][iy][ip];
    }
}
