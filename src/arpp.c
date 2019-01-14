#include "libairs.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Fill data gaps in L2 data. */
void fill_gaps(
  double x[L2_NTRACK][L2_NXTRACK][L2_NLAY],
  double dx,
  double dy);

/* Initialize with AIRS Level-2 data. */
int init_l2(
  airs_l2_t * l2,
  int track,
  int xtrack,
  ctl_t * ctl,
  atm_t * atm);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static atm_t atm, atm_clim;
  static ctl_t ctl;
  static obs_t obs;

  static airs_l1_t l1;
  static airs_l2_t l2;

  FILE *dirlist;

  char cmd[LEN], dirname[LEN];

  double cov = 0, cov_thresh, dx, dy, lat0, lat1, lon0, lon1,
    sza_thresh, z[NP];

  int channel[ND], i, id, ip, iz, nz,
    track, track0, track1, xtrack, xtrack0, xtrack1;

  /* ------------------------------------------------------------
     Read control parameters...
     ------------------------------------------------------------ */

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <airs.nc> <basedir> <dirlist>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read retrieval grid... */
  nz = (int) scan_ctl(argc, argv, "NZ", -1, "", NULL);
  for (iz = 0; iz < nz; iz++)
    z[iz] = scan_ctl(argc, argv, "Z", iz, "", NULL);

  /* Read track range and sampling step... */
  track0 = (int) scan_ctl(argc, argv, "TRACK_MIN", -1, "1", NULL);
  track1 = (int) scan_ctl(argc, argv, "TRACK_MAX", -1, "135", NULL);

  /* Read xtrack range and sampling step... */
  xtrack0 = (int) scan_ctl(argc, argv, "XTRACK_MIN", -1, "1", NULL);
  xtrack1 = (int) scan_ctl(argc, argv, "XTRACK_MAX", -1, "90", NULL);

  /* Read box coordinates... */
  lon0 = scan_ctl(argc, argv, "LON_MIN", -1, "-180", NULL);
  lon1 = scan_ctl(argc, argv, "LON_MAX", -1, "180", NULL);
  lat0 = scan_ctl(argc, argv, "LAT_MIN", -1, "-90", NULL);
  lat1 = scan_ctl(argc, argv, "LAT_MAX", -1, "90", NULL);
  cov_thresh = scan_ctl(argc, argv, "COV_THRESH", -1, "50", NULL);

  /* Smoothing of background... */
  dx = scan_ctl(argc, argv, "DX", -1, "8", NULL);
  dy = scan_ctl(argc, argv, "DY", -1, "2", NULL);

  /* Read SZA threshold... */
  sza_thresh = scan_ctl(argc, argv, "SZA_THRESH", -1, "96", NULL);

  /* ------------------------------------------------------------
     Read AIRS data and initialize...
     ------------------------------------------------------------ */

  /* Read AIRS data... */
  read_l1(argv[2], &l1);
  read_l2(argv[2], &l2);

  /* Check coverage of box... */
  for (track = 0; track < L1_NTRACK; track++)
    for (xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
      if (l1.lon[track][xtrack] >= lon0 && l1.lon[track][xtrack] <= lon1 &&
	  l1.lat[track][xtrack] >= lat0 && l1.lat[track][xtrack] <= lat1)
	cov += 100. / (L1_NTRACK * L1_NXTRACK);
  if (cov < cov_thresh)
    return EXIT_SUCCESS;

  /* Identify radiance channels... */
  for (id = 0; id < ctl.nd; id++) {
    channel[id] = -999;
    for (i = 0; i < L1_NCHAN; i++)
      if (fabs(ctl.nu[id] - l1.nu[i]) < 0.1)
	channel[id] = i;
    if (channel[id] < 0)
      ERRMSG("Cannot identify radiance channel!");
  }

  /* Fill data gaps... */
  fill_gaps(l2.z, dx, dy);
  fill_gaps(l2.t, dx, dy);

  /* Set climatological data for center of granule... */
  atm_clim.np = nz;
  for (iz = 0; iz < nz; iz++) {
    atm_clim.time[iz] = l1.time[L1_NTRACK / 2][L1_NXTRACK / 2];
    atm_clim.z[iz] = z[iz];
    atm_clim.lon[iz] = l1.lon[L1_NTRACK / 2][L1_NXTRACK / 2];
    atm_clim.lat[iz] = l1.lat[L1_NTRACK / 2][L1_NXTRACK / 2];
  }
  climatology(&ctl, &atm_clim);

  /* ------------------------------------------------------------
     Prepare atmospheric data and observation data...
     ------------------------------------------------------------ */

  /* Create directory list... */
  if (!(dirlist = fopen(argv[4], "w")))
    ERRMSG("Cannot create directory list!");

  /* Loop over swaths and scans... */
  for (track = track0 - 1; track < track1; track++)
    for (xtrack = xtrack0 - 1; xtrack < xtrack1; xtrack++) {

      /* Store observation data... */
      obs.nr = 1;
      obs.time[0] = l1.time[track][xtrack];
      obs.obsz[0] = l1.sat_z[track];
      obs.obslon[0] = l1.sat_lon[track];
      obs.obslat[0] = l1.sat_lat[track];
      obs.vpz[0] = 0;
      obs.vplon[0] = l1.lon[track][xtrack];
      obs.vplat[0] = l1.lat[track][xtrack];
      for (id = 0; id < ctl.nd; id++)
	obs.rad[id][0] = l1.rad[track][xtrack][channel[id]];

      /* Flag out 4 micron channels for daytime measurements... */
      if (sza(obs.time[0], obs.obslon[0], obs.obslat[0]) < sza_thresh)
	for (id = 0; id < ctl.nd; id++)
	  if (ctl.nu[id] >= 2000)
	    obs.rad[id][0] = GSL_NAN;

      /* Prepare atmospheric data... */
      copy_atm(&ctl, &atm, &atm_clim, 0);
      for (ip = 0; ip < atm.np; ip++) {
	atm.time[ip] = obs.time[0];
	atm.lon[ip] = obs.vplon[0];
	atm.lat[ip] = obs.vplat[0];
      }

      /* Merge Level-2 data... */
      if (!init_l2(&l2, track, xtrack, &ctl, &atm))
	continue;

      /* Create directory... */
      sprintf(dirname, "%s/swath_%d/scan_%d", argv[3], track + 1, xtrack + 1);
      sprintf(cmd, "mkdir -p %s", dirname);
      if (system(cmd))
	ERRMSG("Cannot create directory!");
      fprintf(dirlist, "%s\n", dirname);

      /* Write observation data... */
      write_obs(dirname, "obs_meas.tab", &ctl, &obs);

      /* Write atmospheric data... */
      write_atm(dirname, "atm_apr.tab", &ctl, &atm);
    }

  /* Close directory list... */
  fclose(dirlist);

  return EXIT_SUCCESS;
}

/************************************************************************/

void fill_gaps(
  double x[L2_NTRACK][L2_NXTRACK][L2_NLAY],
  double dx,
  double dy) {

  static double help[L2_NTRACK][L2_NXTRACK], w, wsum;

  int lay, track, track2, xtrack, xtrack2;

  /* Loop over layers... */
  for (lay = 0; lay < L2_NLAY; lay++) {

    /* Loop over grid points... */
    for (track = 0; track < L2_NTRACK; track++)
      for (xtrack = 0; xtrack < L2_NXTRACK; xtrack++) {

	/* Init... */
	help[track][xtrack] = 0;
	wsum = 0;

	/* Averrage data points... */
	for (track2 = 0; track2 < L2_NTRACK; track2++)
	  for (xtrack2 = 0; xtrack2 < L2_NXTRACK; xtrack2++)
	    if (gsl_finite(x[track2][xtrack2][lay])
		&& x[track2][xtrack2][lay] > 0) {
	      w = exp(-gsl_pow_2((xtrack - xtrack2) / dx)
		      - gsl_pow_2((track - track2) / dy));
	      help[track][xtrack] += w * x[track2][xtrack2][lay];
	      wsum += w;
	    }

	/* Normalize... */
	if (wsum > 0)
	  help[track][xtrack] /= wsum;
	else
	  help[track][xtrack] = GSL_NAN;
      }

    /* Copy grid points... */
    for (track = 0; track < L2_NTRACK; track++)
      for (xtrack = 0; xtrack < L2_NXTRACK; xtrack++)
	x[track][xtrack][lay] = help[track][xtrack];
  }
}

/************************************************************************/

int init_l2(
  airs_l2_t * l2,
  int track,
  int xtrack,
  ctl_t * ctl,
  atm_t * atm) {

  static atm_t atm_airs;

  double k[NW], p, q[NG], t, w, zmax = 0, zmin = 1000;

  int ip, lay;

  /* Reset track- and xtrack-index to match Level-2 data... */
  track /= 3;
  xtrack /= 3;

  /* Store AIRS data in atmospheric data struct... */
  atm_airs.np = 0;
  for (lay = 0; lay < L2_NLAY; lay++)
    if (gsl_finite(l2->z[track][xtrack][lay])) {
      atm_airs.time[atm_airs.np] = l2->time[track][xtrack];
      atm_airs.z[atm_airs.np] = l2->z[track][xtrack][lay];
      atm_airs.lon[atm_airs.np] = l2->lon[track][xtrack];
      atm_airs.lat[atm_airs.np] = l2->lat[track][xtrack];
      atm_airs.p[atm_airs.np] = l2->p[lay];
      atm_airs.t[atm_airs.np] = l2->t[track][xtrack][lay];
      atm_airs.np++;
    }

  /* Check number of levels... */
  if (atm_airs.np <= 0)
    return 0;

  /* Get height range of AIRS data... */
  for (ip = 0; ip < atm_airs.np; ip++) {
    zmax = GSL_MAX(zmax, atm_airs.z[ip]);
    zmin = GSL_MIN(zmin, atm_airs.z[ip]);
  }

  /* Merge AIRS data... */
  for (ip = 0; ip < atm->np; ip++) {

    /* Interpolate AIRS data... */
    intpol_atm(ctl, &atm_airs, atm->z[ip], &p, &t, q, k);

    /* Weighting factor... */
    w = 1;
    if (atm->z[ip] > zmax)
      w = GSL_MAX(1 - (atm->z[ip] - zmax) / 50, 0);
    if (atm->z[ip] < zmin)
      w = GSL_MAX(1 - (zmin - atm->z[ip]) / 50, 0);

    /* Merge... */
    atm->t[ip] = w * t + (1 - w) * atm->t[ip];
    atm->p[ip] = w * p + (1 - w) * atm->p[ip];
  }

  return 1;
}
