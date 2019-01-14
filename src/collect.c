#include "libairs.h"

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/* Number of file types for each variable. */
#define NVAR 7

/* Maximum number of data sets. */
#define NDSMAX 100000

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Add variable attributes to netCDF file. */
void addatt(
  int ncid,
  int varid,
  const char *unit,
  const char *long_name);

/* Create variable in netCDF file. */
void addvar(
  int ncid,
  const char *var_name,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int full,
  int *varid,
  int ndims);

/* Gather retrieval results. */
void analyze(
  ctl_t * ctl,
  const char *dirlistdat,
  size_t * ndir,
  size_t * max_np,
  int *press_full,
  int *temp_full,
  int *q_full,
  int *k_full);

/* Read atmospheric data. */
void collect_read_atm(
  char *dirname,
  char *filename,
  ctl_t * ctl,
  atm_t * atm_dest,
  atm_t * atm_src);

/* Write variable data to netcdf file. */
void writevar(
  int ncid,
  int varid[],
  int full,
  double *atm,
  double *tot,
  double *noi,
  double *fm,
  double *cont,
  double *res,
  double *apr,
  size_t ind_start[],
  size_t ind_end[]);

/* Read cost function data. */
double read_chisq(
  const char *dirname);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static atm_t atm, atm_noise, atm_fm, atm_tot, atm_cont, atm_res, atm_apr;
  static ctl_t ctl;

  FILE *dirlist;

  double chisq;

  char dir[LEN], longname[LEN], varname[LEN];

  int chisq_id[NVAR], cnt = -1, dimid[2], dimid_chisq[1], ig, iw,
    k_full[ND] = { 0 }, k_id[ND][NVAR], lat_id[NVAR], lon_id[NVAR],
    na[NDSMAX], ncid, nds_id, np_id, np_dimid, nds_dimid,
    press_full = 0, press_id[NVAR], q_full[NG] = {
  0}, q_id[NG][NVAR], temp_full = 0, temp_id[NVAR], time_id[NVAR], z_id[NVAR];

  size_t i, ind_start[] = { 0, 0 }, ind_end[] = {
  0, 0}, ind_start_chisq[] = {
  0}, ind_end_chisq[] = {
  0}, max_np, ndir = 0;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <dirlist> <out.nc>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Get maximum size of variable arrays... */
  analyze(&ctl, argv[2], &ndir, &max_np,
	  &press_full, &temp_full, q_full, k_full);

  /* Check number of data sets... */
  if (ndir >= NDSMAX)
    ERRMSG("Too many data sets!");

  /* Create netCDF file... */
  NC(nc_create(argv[3], NC_CLOBBER, &ncid));

  /* Set dimensions... */
  NC(nc_def_dim(ncid, "np", max_np, &np_dimid));
  NC(nc_def_dim(ncid, "nds", ndir, &nds_dimid));

  /* Set dimids... */
  dimid[0] = np_dimid;
  dimid[1] = nds_dimid;
  dimid_chisq[0] = nds_dimid;

  /* Coordinate variables... */
  addvar(ncid, "nds", "1", "data set index",
	 NC_INT, dimid_chisq, 0, &nds_id, 1);
  addvar(ncid, "np", "1", "data point index", NC_INT, dimid, 0, &np_id, 1);

  /* Atmospheric data... */
  addvar(ncid, "time", "s", "time (seconds since 2000-01-01T00:00Z)",
	 NC_DOUBLE, dimid, 0, time_id, 2);
  addvar(ncid, "z", "km", "altitude", NC_DOUBLE, dimid, 0, z_id, 2);
  addvar(ncid, "lon", "deg", "longitude", NC_DOUBLE, dimid, 0, lon_id, 2);
  addvar(ncid, "lat", "deg", "latitude", NC_DOUBLE, dimid, 0, lat_id, 2);
  addvar(ncid, "press", "hPa", "pressure",
	 NC_FLOAT, dimid, press_full, press_id, 2);
  addvar(ncid, "temp", "K", "temperature",
	 NC_FLOAT, dimid, temp_full, temp_id, 2);
  for (ig = 0; ig < ctl.ng; ig++) {
    sprintf(longname, "%s %s", ctl.emitter[ig], "volume mixing ratio");
    addvar(ncid, ctl.emitter[ig], "1", longname,
	   NC_FLOAT, dimid, q_full[ig], q_id[ig], 2);
  }
  for (iw = 0; iw < ctl.nw; iw++) {
    sprintf(varname, "extinct_win%d", iw);
    sprintf(longname, "extinction (window %d)", iw);
    addvar(ncid, varname, "km^-1", longname,
	   NC_FLOAT, dimid, k_full[iw], k_id[iw], 2);
  }

  /* Cost function data... */
  addvar(ncid, "chisq", "1", "normalized cost function",
	 NC_FLOAT, dimid_chisq, 0, chisq_id, 1);

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Write coordinate variables... */
  ind_start[0] = 0;
  ind_end[0] = ndir;
  for (i = 0; i < ind_end[0]; i++)
    na[i] = (int) i + 1;
  NC(nc_put_vara_int(ncid, nds_id, ind_start, ind_end, na));

  ind_start[0] = 0;
  ind_end[0] = max_np;
  for (i = 0; i < ind_end[0]; i++)
    na[i] = (int) i + 1;
  NC(nc_put_vara_int(ncid, np_id, ind_start, ind_end, na));

  /* Open directory list... */
  if (!(dirlist = fopen(argv[2], "r")))
    ERRMSG("Cannot open directory list!");

  /* Loop over directories... */
  while (fscanf(dirlist, "%s", dir) != EOF) {

    /* Increment counter... */
    cnt++;

    /* Read atmospheric data... */
    read_atm(dir, "atm_apr.tab", &ctl, &atm_apr);
    collect_read_atm(dir, "atm_final.tab", &ctl, &atm, &atm_apr);
    collect_read_atm(dir, "atm_err_total.tab", &ctl, &atm_tot, &atm_apr);
    collect_read_atm(dir, "atm_err_noise.tab", &ctl, &atm_noise, &atm_apr);
    collect_read_atm(dir, "atm_err_formod.tab", &ctl, &atm_fm, &atm_apr);
    collect_read_atm(dir, "atm_cont.tab", &ctl, &atm_cont, &atm_apr);
    collect_read_atm(dir, "atm_res.tab", &ctl, &atm_res, &atm_apr);
    chisq = read_chisq(dir);

    /* Set indices... */
    ind_start[0] = 0;
    ind_start[1] = (size_t) cnt;
    ind_end[0] = (size_t) atm.np;
    ind_end[1] = 1;
    ind_start_chisq[0] = (size_t) cnt;
    ind_end_chisq[0] = 1;

    /* Write atmospheric data... */
    writevar(ncid, time_id, 0, atm.time, NULL, NULL, NULL, NULL, NULL, NULL,
	     ind_start, ind_end);
    writevar(ncid, z_id, 0, atm.z, NULL, NULL, NULL, NULL, NULL, NULL,
	     ind_start, ind_end);
    writevar(ncid, lon_id, 0, atm.lon, NULL, NULL, NULL, NULL, NULL, NULL,
	     ind_start, ind_end);
    writevar(ncid, lat_id, 0, atm.lat, NULL, NULL, NULL, NULL, NULL, NULL,
	     ind_start, ind_end);
    writevar(ncid, press_id, press_full, atm.p, atm_tot.p, atm_noise.p,
	     atm_fm.p, atm_cont.p, atm_res.p, atm_apr.p, ind_start, ind_end);
    writevar(ncid, temp_id, temp_full, atm.t, atm_tot.t, atm_noise.t,
	     atm_fm.t, atm_cont.t, atm_res.t, atm_apr.t, ind_start, ind_end);
    for (ig = 0; ig < ctl.ng; ig++)
      writevar(ncid, q_id[ig], q_full[ig], atm.q[ig], atm_tot.q[ig],
	       atm_noise.q[ig], atm_fm.q[ig], atm_cont.q[ig],
	       atm_res.q[ig], atm_apr.q[ig], ind_start, ind_end);
    for (iw = 0; iw < ctl.nw; iw++)
      writevar(ncid, k_id[iw], k_full[iw], atm.k[iw], atm_tot.k[iw],
	       atm_noise.k[iw], atm_fm.k[iw], atm_cont.k[iw], atm_res.k[iw],
	       atm_apr.k[iw], ind_start, ind_end);

    /* Write cost function data... */
    writevar(ncid, chisq_id, 0, &chisq, NULL, NULL, NULL, NULL, NULL, NULL,
	     ind_start_chisq, ind_end_chisq);
  }

  /* Close file... */
  NC(nc_close(ncid));

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void addatt(
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

void addvar(
  int ncid,
  const char *var_name,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int full,
  int *varid,
  int ndims) {

  char lname[LEN], var[LEN];

  /* Retrieval results... */
  NC(nc_def_var(ncid, var_name, type, ndims, dimid, &varid[0]));
  sprintf(lname, "%s", longname);
  addatt(ncid, varid[0], unit, lname);

  /* Write full output? */
  if (full) {

    /* Total error... */
    sprintf(var, "%s_total", var_name);
    NC(nc_def_var(ncid, var, NC_FLOAT, 2, dimid, &varid[1]));
    sprintf(lname, "%s total error", longname);
    addatt(ncid, varid[1], unit, lname);

    /* Noise error... */
    sprintf(var, "%s_noise", var_name);
    NC(nc_def_var(ncid, var, NC_FLOAT, 2, dimid, &varid[2]));
    sprintf(lname, "%s noise error", longname);
    addatt(ncid, varid[2], unit, lname);

    /* Forward model error... */
    sprintf(var, "%s_formod", var_name);
    NC(nc_def_var(ncid, var, NC_FLOAT, 2, dimid, &varid[3]));
    sprintf(lname, "%s forward model error", longname);
    addatt(ncid, varid[3], unit, lname);

    /* Measurement content... */
    sprintf(var, "%s_cont", var_name);
    NC(nc_def_var(ncid, var, NC_FLOAT, 2, dimid, &varid[4]));
    sprintf(lname, "%s measurement content", longname);
    addatt(ncid, varid[4], "1", lname);

    /* Resolution... */
    sprintf(var, "%s_res", var_name);
    NC(nc_def_var(ncid, var, NC_FLOAT, 2, dimid, &varid[5]));
    sprintf(lname, "%s resolution", longname);
    addatt(ncid, varid[5], "1", lname);

    /* A priori data... */
    sprintf(var, "%s_apr", var_name);
    NC(nc_def_var(ncid, var, NC_FLOAT, 2, dimid, &varid[6]));
    sprintf(lname, "%s a priori", longname);
    addatt(ncid, varid[6], unit, lname);
  }
}

/*****************************************************************************/

void analyze(
  ctl_t * ctl,
  const char *dirlistdat,
  size_t * ndir,
  size_t * max_np,
  int *press_full,
  int *temp_full,
  int *q_full,
  int *k_full) {

  static atm_t atm;

  FILE *dirlist;

  char dir[LEN];

  int ig, ip, iw;

  /* Initialize... */
  *ndir = *max_np = 0;
  *press_full = *temp_full = 0;
  for (ig = 0; ig < ctl->ng; ig++)
    q_full[ig] = 0;
  for (iw = 0; iw < ctl->nw; iw++)
    k_full[iw] = 0;

  /* Open directory list... */
  if (!(dirlist = fopen(dirlistdat, "r")))
    ERRMSG("Cannot open directory list!");

  /* Loop over directories... */
  while (fscanf(dirlist, "%s", dir) != EOF) {

    /* Increment profile counter... */
    ++(*ndir);

    /* Read atmospheric data... */
    read_atm(dir, "atm_apr.tab", ctl, &atm);

    /* Analyze atmospheric data... */
    *max_np = GSL_MAX((size_t) atm.np, *max_np);
    for (ip = 0; ip < atm.np; ip++) {
      *press_full = *press_full ||
	(atm.z[ip] >= ctl->retp_zmin && atm.z[ip] <= ctl->retp_zmax);
      *temp_full = *temp_full ||
	(atm.z[ip] >= ctl->rett_zmin && atm.z[ip] <= ctl->rett_zmax);
      for (ig = 0; ig < ctl->ng; ig++)
	q_full[ig] = q_full[ig] || (atm.z[ip] >= ctl->retq_zmin[ig]
				    && atm.z[ip] <= ctl->retq_zmax[ig]);
      for (iw = 0; iw < ctl->nw; iw++)
	k_full[iw] = k_full[iw] || (atm.z[ip] >= ctl->retk_zmin[iw]
				    && atm.z[ip] <= ctl->retk_zmax[iw]);
    }
  }

  /* Close directory list... */
  fclose(dirlist);
}

/*****************************************************************************/

void collect_read_atm(
  char *dirname,
  char *filename,
  ctl_t * ctl,
  atm_t * atm_dest,
  atm_t * atm_src) {

  FILE *in;

  char file[LEN];

  /* Set filename... */
  sprintf(file, "%s/%s", dirname, filename);

  /* Try to read file... */
  if (!(in = fopen(file, "r")))
    copy_atm(ctl, atm_dest, atm_src, 1);
  else {
    read_atm(dirname, filename, ctl, atm_dest);
    fclose(in);
  }
}

/*****************************************************************************/

void writevar(
  int ncid,
  int varid[],
  int full,
  double *atm,
  double *tot,
  double *noi,
  double *fm,
  double *cont,
  double *res,
  double *apr,
  size_t ind_start[],
  size_t ind_end[]) {

  /* Retrieval results... */
  NC(nc_put_vara_double(ncid, varid[0], ind_start, ind_end, atm));

  /* Write full output? */
  if (full) {

    /* Total error... */
    NC(nc_put_vara_double(ncid, varid[1], ind_start, ind_end, tot));

    /* Noise error... */
    NC(nc_put_vara_double(ncid, varid[2], ind_start, ind_end, noi));

    /* Forward model error... */
    NC(nc_put_vara_double(ncid, varid[3], ind_start, ind_end, fm));

    /* Measurement content... */
    NC(nc_put_vara_double(ncid, varid[4], ind_start, ind_end, cont));

    /* Resolution... */
    NC(nc_put_vara_double(ncid, varid[5], ind_start, ind_end, res));

    /* A priori... */
    NC(nc_put_vara_double(ncid, varid[6], ind_start, ind_end, apr));
  }
}

/*****************************************************************************/

double read_chisq(
  const char *dirname) {

  FILE *in;

  char file[LEN], line[LEN], *tok;

  double chisq, dummy;

  /* Set filename... */
  sprintf(file, "%s/costs.tab", dirname);

  /* Write info... */
  printf("Read cost function data: %s\n", file);

  /* Open file... */
  if (!(in = fopen(file, "r")))
    return GSL_NAN;

  /* Read data... */
  while (fgets(line, LEN, in)) {
    TOK(line, tok, "%lg", dummy);
    TOK(NULL, tok, "%lg", chisq);
  }

  /* Close file... */
  fclose(in);

  return chisq;
}
