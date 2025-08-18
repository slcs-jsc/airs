/*
  This file is part of the AIRS Code Collection.
  
  the AIRS Code Collections is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.
  
  The AIRS Code Collection is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with the AIRS Code Collection. If not, see
  <http://www.gnu.org/licenses/>.
  
  Copyright (C) 2019-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Retrieval processor for AIRS.
*/

#include <mpi.h>
#include <omp.h>
#include <netcdf.h>
#include "jurassic.h"

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Execute netCDF library command and check result. */
#define NC(cmd) {				     \
  int nc_result=(cmd);				     \
  if(nc_result!=NC_NOERR)			     \
    ERRMSG("%s", nc_strerror(nc_result));	     \
}

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Number of AIRS radiance channels (don't change). */
#define L1_NCHAN 34

/*! Along-track size of AIRS radiance granule (don't change). */
#define L1_NTRACK 135

/*! Across-track size of AIRS radiance granule (don't change). */
#define L1_NXTRACK 90

/*! Number of AIRS pressure layers (don't change). */
#define L2_NLAY 27

/*! Along-track size of AIRS retrieval granule (don't change). */
#define L2_NTRACK 45

/*! Across-track size of AIRS retrieval granule (don't change). */
#define L2_NXTRACK 30

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! Buffer for netCDF data. */
typedef struct {

  /*! NetCDF file ID. */
  int ncid;

  /*! Number of retrieval altitudes. */
  int np;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double l1_time[L1_NTRACK][L1_NXTRACK];

  /*! Footprint longitude [deg]. */
  double l1_lon[L1_NTRACK][L1_NXTRACK];

  /*! Footprint latitude [deg]. */
  double l1_lat[L1_NTRACK][L1_NXTRACK];

  /*! Satellite altitude [km]. */
  double l1_sat_z[L1_NTRACK];

  /*! Satellite longitude [deg]. */
  double l1_sat_lon[L1_NTRACK];

  /*! Satellite latitude [deg]. */
  double l1_sat_lat[L1_NTRACK];

  /*! Channel frequencies [cm^-1]. */
  double l1_nu[L1_NCHAN];

  /*! Radiance [W/(m^2 sr cm^-1)]. */
  float l1_rad[L1_NTRACK][L1_NXTRACK][L1_NCHAN];

  /*! Altitude [km]. */
  double l2_z[L2_NTRACK][L2_NXTRACK][L2_NLAY];

  /*! Pressure [hPa]. */
  double l2_p[L2_NLAY];

  /*! Temperature [K]. */
  double l2_t[L2_NTRACK][L2_NXTRACK][L2_NLAY];

  /*! Altitude [km]. */
  float ret_z[NP];

  /*! Pressure [hPa]. */
  float ret_p[L1_NTRACK * L1_NXTRACK];

  /*! Temperature [K]. */
  float ret_t[L1_NTRACK * L1_NXTRACK * NP];

} ncd_t;

/*! Retrieval control parameters. */
typedef struct {

  /*! Recomputation of kernel matrix (number of iterations). */
  int kernel_recomp;

  /*! Maximum number of iterations. */
  int conv_itmax;

  /*! Minimum normalized step size in state space. */
  double conv_dmin;

  /*! Forward model error [%]. */
  double err_formod[ND];

  /*! Noise error [W/(m^2 sr cm^-1)]. */
  double err_noise[ND];

  /*! Pressure error [%]. */
  double err_press;

  /*! Vertical correlation length for pressure error [km]. */
  double err_press_cz;

  /*! Horizontal correlation length for pressure error [km]. */
  double err_press_ch;

  /*! Temperature error [K]. */
  double err_temp;

  /*! Vertical correlation length for temperature error [km]. */
  double err_temp_cz;

  /*! Horizontal correlation length for temperature error [km]. */
  double err_temp_ch;

  /*! Volume mixing ratio error [%]. */
  double err_q[NG];

  /*! Vertical correlation length for volume mixing ratio error [km]. */
  double err_q_cz[NG];

  /*! Horizontal correlation length for volume mixing ratio error [km]. */
  double err_q_ch[NG];

  /*! Extinction error [1/km]. */
  double err_k[NW];

  /*! Vertical correlation length for extinction error [km]. */
  double err_k_cz[NW];

  /*! Horizontal correlation length for extinction error [km]. */
  double err_k_ch[NW];

} ret_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Create variable in netCDF file. */
void add_var(
  int ncid,
  const char *varname,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int *varid,
  int ndims);

/*! Buffer netCDF data. */
void buffer_nc(
  atm_t * atm,
  double chisq,
  ncd_t * ncd,
  int track,
  int xtrack,
  int np0,
  int np1);

/*! Compute cost function. */
double cost_function(
  gsl_vector * dx,
  gsl_vector * dy,
  gsl_matrix * s_a_inv,
  gsl_vector * sig_eps_inv);

/*! Fill data gaps in L2 data. */
void fill_gaps(
  double x[L2_NTRACK][L2_NXTRACK][L2_NLAY],
  double cx,
  double cy);

/*! Initialize with AIRS Level-2 data. */
void init_l2(
  ncd_t * ncd,
  int track,
  int xtrack,
  ctl_t * ctl,
  atm_t * atm);

/*! Invert symmetric matrix. */
void matrix_invert(
  gsl_matrix * a);

/*! Compute matrix product A^TBA or ABA^T for diagonal matrix B. */
void matrix_product(
  gsl_matrix * a,
  gsl_vector * b,
  int transpose,
  gsl_matrix * c);

/*! Carry out optimal estimation retrieval. */
void optimal_estimation(
  ret_t * ret,
  ctl_t * ctl,
  tbl_t * tbl,
  obs_t * obs_meas,
  obs_t * obs_i,
  atm_t * atm_apr,
  atm_t * atm_i,
  double *chisq);

/*! Read netCDF file. */
void read_nc(
  char *filename,
  ncd_t * ncd);

/*! Read retrieval control parameters. */
void read_ret_ctl(
  int argc,
  char *argv[],
  ctl_t * ctl,
  ret_t * ret);

/*! Set a priori covariance. */
void set_cov_apr(
  ret_t * ret,
  ctl_t * ctl,
  atm_t * atm,
  int *iqa,
  int *ipa,
  gsl_matrix * s_a);

/*! Set measurement errors. */
void set_cov_meas(
  ret_t * ret,
  ctl_t * ctl,
  obs_t * obs,
  gsl_vector * sig_noise,
  gsl_vector * sig_formod,
  gsl_vector * sig_eps_inv);

/*! Write to netCDF file... */
void write_nc(
  char *filename,
  ncd_t * ncd);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static atm_t atm_apr, atm_clim, atm_i;
  static obs_t obs_i, obs_meas;
  static ncd_t ncd;
  static ret_t ret;

  FILE *in;

  char filename[LEN];

  double chisq, chisq_min, chisq_max, chisq_mean, z[NP];

  int channel[ND], m, ntask = -1, rank, size;

  /* ------------------------------------------------------------
     Init...
     ------------------------------------------------------------ */

  /* MPI... */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Measure CPU time... */
  TIMER("total", 1);

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <filelist>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  read_ret_ctl(argc, argv, &ctl, &ret);

  /* Initialize look-up tables... */
  tbl_t *tbl = read_tbl(&ctl);

  /* Read retrieval grid... */
  const int nz = (int) scan_ctl(argc, argv, "NZ", -1, "", NULL);
  if (nz > NP)
    ERRMSG("Too many altitudes!");
  for (int iz = 0; iz < nz; iz++)
    z[iz] = scan_ctl(argc, argv, "Z", iz, "", NULL);

  /* Read track range... */
  const int track0 = (int) scan_ctl(argc, argv, "TRACK_MIN", -1, "0", NULL);
  const int track1 = (int) scan_ctl(argc, argv, "TRACK_MAX", -1, "134", NULL);

  /* Read xtrack range... */
  const int xtrack0 = (int) scan_ctl(argc, argv, "XTRACK_MIN", -1, "0", NULL);
  const int xtrack1 =
    (int) scan_ctl(argc, argv, "XTRACK_MAX", -1, "89", NULL);

  /* Read height range... */
  int np0 = (int) scan_ctl(argc, argv, "NP_MIN", -1, "0", NULL);
  int np1 = (int) scan_ctl(argc, argv, "NP_MAX", -1, "100", NULL);
  np1 = GSL_MIN(np1, nz - 1);

  /* Background smoothing... */
  const double sx = scan_ctl(argc, argv, "SX", -1, "8", NULL);
  const double sy = scan_ctl(argc, argv, "SY", -1, "2", NULL);

  /* SZA threshold... */
  const double sza_thresh = scan_ctl(argc, argv, "SZA", -1, "96", NULL);

  /* ------------------------------------------------------------
     Distribute granules...
     ------------------------------------------------------------ */

  /* Open filelist... */
  printf("Read filelist: %s\n", argv[2]);
  if (!(in = fopen(argv[2], "r")))
    ERRMSG("Cannot open filelist!");

  /* Loop over netCDF files... */
  while (fscanf(in, "%s", filename) != EOF) {

    /* Distribute files with MPI... */
    if ((++ntask) % size != rank)
      continue;

    /* Write info... */
    printf("Retrieve file %s on rank %d of %d (with %d threads)...\n",
	   filename, rank + 1, size, omp_get_max_threads());

    /* ------------------------------------------------------------
       Initialize retrieval...
       ------------------------------------------------------------ */

    /* Read netCDF file... */
    read_nc(filename, &ncd);

    /* Identify radiance channels... */
    for (int id = 0; id < ctl.nd; id++) {
      channel[id] = -999;
      for (int i = 0; i < L1_NCHAN; i++)
	if (fabs(ctl.nu[id] - ncd.l1_nu[i]) < 0.1)
	  channel[id] = i;
      if (channel[id] < 0)
	ERRMSG("Cannot identify radiance channel!");
    }

    /* Fill data gaps... */
    fill_gaps(ncd.l2_t, sx, sy);
    fill_gaps(ncd.l2_z, sx, sy);

    /* Set climatological data for center of granule... */
    atm_clim.np = nz;
    for (int iz = 0; iz < nz; iz++)
      atm_clim.z[iz] = z[iz];
    climatology(&ctl, &atm_clim);

    /* ------------------------------------------------------------
       Retrieval...
       ------------------------------------------------------------ */

    /* Get chi^2 statistics... */
    chisq_min = 1e100;
    chisq_max = -1e100;
    chisq_mean = 0;
    m = 0;

    /* Loop over swaths... */
    for (int track = track0; track <= track1; track++) {

      /* Measure CPU time... */
      TIMER("retrieval", 1);

      /* Loop over scan... */
      for (int xtrack = xtrack0; xtrack <= xtrack1; xtrack++) {

	/* Store observation data... */
	obs_meas.nr = 1;
	obs_meas.time[0] = ncd.l1_time[track][xtrack];
	obs_meas.obsz[0] = ncd.l1_sat_z[track];
	obs_meas.obslon[0] = ncd.l1_sat_lon[track];
	obs_meas.obslat[0] = ncd.l1_sat_lat[track];
	obs_meas.vplon[0] = ncd.l1_lon[track][xtrack];
	obs_meas.vplat[0] = ncd.l1_lat[track][xtrack];
	for (int id = 0; id < ctl.nd; id++)
	  obs_meas.rad[id][0] = ncd.l1_rad[track][xtrack][channel[id]];

	/* Flag out 4 micron channels for daytime measurements... */
	if (sza(obs_meas.time[0], obs_meas.obslon[0], obs_meas.obslat[0])
	    < sza_thresh)
	  for (int id = 0; id < ctl.nd; id++)
	    if (ctl.nu[id] >= 2000)
	      obs_meas.rad[id][0] = GSL_NAN;

	/* Prepare atmospheric data... */
	copy_atm(&ctl, &atm_apr, &atm_clim, 0);
	for (int ip = 0; ip < atm_apr.np; ip++) {
	  atm_apr.time[ip] = obs_meas.time[0];
	  atm_apr.lon[ip] = obs_meas.vplon[0];
	  atm_apr.lat[ip] = obs_meas.vplat[0];
	}

	/* Merge Level-2 data... */
	init_l2(&ncd, track, xtrack, &ctl, &atm_apr);

	/* Retrieval... */
	optimal_estimation(&ret, &ctl, tbl, &obs_meas, &obs_i,
			   &atm_apr, &atm_i, &chisq);

	/* Get chi^2 statistics... */
	if (gsl_finite(chisq)) {
	  chisq_min = GSL_MIN(chisq_min, chisq);
	  chisq_max = GSL_MAX(chisq_max, chisq);
	  chisq_mean += chisq;
	  m++;
	}

	/* Buffer results... */
	buffer_nc(&atm_i, chisq, &ncd, track, xtrack, np0, np1);
      }

      /* Measure CPU time... */
      TIMER("retrieval", 3);
    }

    /* ------------------------------------------------------------
       Finalize...
       ------------------------------------------------------------ */

    /* Write netCDF file... */
    write_nc(filename, &ncd);

    /* Write info... */
    printf("chi^2: min= %g / mean= %g / max= %g / m= %d\n",
	   chisq_min, chisq_mean / m, chisq_max, m);
    printf("Retrieval finished on rank %d of %d!\n", rank, size);
  }

  /* Close file list... */
  fclose(in);

  /* Measure CPU time... */
  TIMER("total", 3);

  /* Report memory usage... */
  printf("MEMORY_ATM = %g MByte\n", 4. * sizeof(atm_t) / 1024. / 1024.);
  printf("MEMORY_CTL = %g MByte\n", 1. * sizeof(ctl_t) / 1024. / 1024.);
  printf("MEMORY_NCD = %g MByte\n", 1. * sizeof(ncd_t) / 1024. / 1024.);
  printf("MEMORY_OBS = %g MByte\n", 3. * sizeof(atm_t) / 1024. / 1024.);
  printf("MEMORY_RET = %g MByte\n", 1. * sizeof(ret_t) / 1024. / 1024.);
  printf("MEMORY_TBL = %g MByte\n", 1. * sizeof(tbl_t) / 1024. / 1024.);

  /* Report problem size... */
  printf("SIZE_TASKS = %d\n", size);
  printf("SIZE_THREADS = %d\n", omp_get_max_threads());

  /* MPI... */
  MPI_Finalize();

  return EXIT_SUCCESS;
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

void buffer_nc(
  atm_t *atm,
  double chisq,
  ncd_t *ncd,
  int track,
  int xtrack,
  int np0,
  int np1) {

  /* Set number of data points... */
  ncd->np = np1 - np0 + 1;

  /* Save retrieval data... */
  for (int ip = np0; ip <= np1; ip++) {
    ncd->ret_z[ip - np0] = (float) atm->z[ip];
    ncd->ret_p[track * L1_NXTRACK + xtrack] = (float) atm->p[np0];
    ncd->ret_t[(track * L1_NXTRACK + xtrack) * ncd->np + ip - np0] =
      (gsl_finite(chisq) ? (float) atm->t[ip] : GSL_NAN);
  }
}

/*****************************************************************************/

double cost_function(
  gsl_vector *dx,
  gsl_vector *dy,
  gsl_matrix *s_a_inv,
  gsl_vector *sig_eps_inv) {

  double chisq_a, chisq_m = 0;

  /* Get sizes... */
  const size_t m = dy->size;
  const size_t n = dx->size;

  /* Allocate... */
  gsl_vector *x_aux = gsl_vector_alloc(n);
  gsl_vector *y_aux = gsl_vector_alloc(m);

  /* Determine normalized cost function...
     (chi^2 = 1/m * [dy^T * S_eps^{-1} * dy + dx^T * S_a^{-1} * dx]) */
  for (size_t i = 0; i < m; i++)
    chisq_m +=
      gsl_pow_2(gsl_vector_get(dy, i) * gsl_vector_get(sig_eps_inv, i));
  gsl_blas_dgemv(CblasNoTrans, 1.0, s_a_inv, dx, 0.0, x_aux);
  gsl_blas_ddot(dx, x_aux, &chisq_a);

  /* Free... */
  gsl_vector_free(x_aux);
  gsl_vector_free(y_aux);

  /* Return cost function value... */
  return (chisq_m + chisq_a) / (double) m;
}

/************************************************************************/

void fill_gaps(
  double x[L2_NTRACK][L2_NXTRACK][L2_NLAY],
  double cx,
  double cy) {

  double help[L2_NTRACK][L2_NXTRACK];

  /* Loop over layers... */
  for (int lay = 0; lay < L2_NLAY; lay++) {

    /* Loop over grid points... */
    for (int track = 0; track < L2_NTRACK; track++)
      for (int xtrack = 0; xtrack < L2_NXTRACK; xtrack++) {

	/* Init... */
	help[track][xtrack] = 0;
	double wsum = 0;

	/* Averrage data points... */
	for (int track2 = 0; track2 < L2_NTRACK; track2++)
	  for (int xtrack2 = 0; xtrack2 < L2_NXTRACK; xtrack2++)
	    if (gsl_finite(x[track2][xtrack2][lay])
		&& x[track2][xtrack2][lay] > 0) {
	      const double w = exp(-gsl_pow_2((xtrack - xtrack2) / cx)
				   - gsl_pow_2((track - track2) / cy));
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
    for (int track = 0; track < L2_NTRACK; track++)
      for (int xtrack = 0; xtrack < L2_NXTRACK; xtrack++)
	x[track][xtrack][lay] = help[track][xtrack];
  }
}

/************************************************************************/

void init_l2(
  ncd_t *ncd,
  int track,
  int xtrack,
  ctl_t *ctl,
  atm_t *atm) {

  static atm_t atm_airs;

  double k[NW], p, q[NG], t, zmax = 0, zmin = 1000;

  /* Reset track- and xtrack-index to match Level-2 data... */
  track /= 3;
  xtrack /= 3;

  /* Store AIRS data in atmospheric data struct... */
  atm_airs.np = 0;
  for (int lay = 0; lay < L2_NLAY; lay++)
    if (gsl_finite(ncd->l2_z[track][xtrack][lay])) {
      atm_airs.z[atm_airs.np] = ncd->l2_z[track][xtrack][lay];
      atm_airs.p[atm_airs.np] = ncd->l2_p[lay];
      atm_airs.t[atm_airs.np] = ncd->l2_t[track][xtrack][lay];
      if ((++atm_airs.np) > NP)
	ERRMSG("Too many layers!");
    }

  /* Check number of levels... */
  if (atm_airs.np <= 0)
    return;

  /* Get height range of AIRS data... */
  for (int ip = 0; ip < atm_airs.np; ip++) {
    zmax = GSL_MAX(zmax, atm_airs.z[ip]);
    zmin = GSL_MIN(zmin, atm_airs.z[ip]);
  }

  /* Merge AIRS data... */
  for (int ip = 0; ip < atm->np; ip++) {

    /* Interpolate AIRS data... */
    intpol_atm(ctl, &atm_airs, atm->z[ip], &p, &t, q, k);

    /* Weighting factor... */
    double w = 1;
    if (atm->z[ip] > zmax)
      w = GSL_MAX(1 - (atm->z[ip] - zmax) / 50, 0);
    if (atm->z[ip] < zmin)
      w = GSL_MAX(1 - (zmin - atm->z[ip]) / 50, 0);

    /* Merge... */
    atm->t[ip] = w * t + (1 - w) * atm->t[ip];
    atm->p[ip] = w * p + (1 - w) * atm->p[ip];
  }
}

/*****************************************************************************/

void matrix_invert(
  gsl_matrix *a) {

  size_t diag = 1;

  /* Get size... */
  const size_t n = a->size1;

  /* Check if matrix is diagonal... */
  for (size_t i = 0; i < n && diag; i++)
    for (size_t j = i + 1; j < n; j++)
      if (gsl_matrix_get(a, i, j) != 0) {
	diag = 0;
	break;
      }

  /* Quick inversion of diagonal matrix... */
  if (diag)
    for (size_t i = 0; i < n; i++)
      gsl_matrix_set(a, i, i, 1 / gsl_matrix_get(a, i, i));

  /* Matrix inversion by means of Cholesky decomposition... */
  else {
    gsl_linalg_cholesky_decomp(a);
    gsl_linalg_cholesky_invert(a);
  }
}

/*****************************************************************************/

void matrix_product(
  gsl_matrix *a,
  gsl_vector *b,
  int transpose,
  gsl_matrix *c) {

  /* Set sizes... */
  const size_t m = a->size1;
  const size_t n = a->size2;

  /* Allocate... */
  gsl_matrix *aux = gsl_matrix_alloc(m, n);

  /* Compute A^T B A... */
  if (transpose == 1) {

    /* Compute B^1/2 A... */
    for (size_t i = 0; i < m; i++)
      for (size_t j = 0; j < n; j++)
	gsl_matrix_set(aux, i, j,
		       gsl_vector_get(b, i) * gsl_matrix_get(a, i, j));

    /* Compute A^T B A = (B^1/2 A)^T (B^1/2 A)... */
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, aux, aux, 0.0, c);
  }

  /* Compute A B A^T... */
  else if (transpose == 2) {

    /* Compute A B^1/2... */
    for (size_t i = 0; i < m; i++)
      for (size_t j = 0; j < n; j++)
	gsl_matrix_set(aux, i, j,
		       gsl_matrix_get(a, i, j) * gsl_vector_get(b, j));

    /* Compute A B A^T = (A B^1/2) (A B^1/2)^T... */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, aux, aux, 0.0, c);
  }

  /* Free... */
  gsl_matrix_free(aux);
}

/*****************************************************************************/

void optimal_estimation(
  ret_t *ret,
  ctl_t *ctl,
  tbl_t *tbl,
  obs_t *obs_meas,
  obs_t *obs_i,
  atm_t *atm_apr,
  atm_t *atm_i,
  double *chisq) {

  static int ipa[N], iqa[N];

  double disq = 0, lmpar = 0.001;

  /* ------------------------------------------------------------
     Initialize...
     ------------------------------------------------------------ */

  /* Get sizes... */
  const size_t m = obs2y(ctl, obs_meas, NULL, NULL, NULL);
  const size_t n = atm2x(ctl, atm_apr, NULL, iqa, ipa);
  if (m == 0 || n == 0) {
    *chisq = GSL_NAN;
    return;
  }

  /* Allocate... */
  gsl_matrix *a = gsl_matrix_alloc(n, n);
  gsl_matrix *cov = gsl_matrix_alloc(n, n);
  gsl_matrix *k_i = gsl_matrix_alloc(m, n);
  gsl_matrix *s_a_inv = gsl_matrix_alloc(n, n);

  gsl_vector *b = gsl_vector_alloc(n);
  gsl_vector *dx = gsl_vector_alloc(n);
  gsl_vector *dy = gsl_vector_alloc(m);
  gsl_vector *sig_eps_inv = gsl_vector_alloc(m);
  gsl_vector *sig_formod = gsl_vector_alloc(m);
  gsl_vector *sig_noise = gsl_vector_alloc(m);
  gsl_vector *x_a = gsl_vector_alloc(n);
  gsl_vector *x_i = gsl_vector_alloc(n);
  gsl_vector *x_step = gsl_vector_alloc(n);
  gsl_vector *y_aux = gsl_vector_alloc(m);
  gsl_vector *y_i = gsl_vector_alloc(m);
  gsl_vector *y_m = gsl_vector_alloc(m);

  /* Set initial state... */
  copy_atm(ctl, atm_i, atm_apr, 0);
  copy_obs(ctl, obs_i, obs_meas, 0);
  formod(ctl, tbl, atm_i, obs_i);

  /* Set state vectors and observation vectors... */
  atm2x(ctl, atm_apr, x_a, NULL, NULL);
  atm2x(ctl, atm_i, x_i, NULL, NULL);
  obs2y(ctl, obs_meas, y_m, NULL, NULL);
  obs2y(ctl, obs_i, y_i, NULL, NULL);

  /* Set inverse a priori covariance S_a^-1... */
  set_cov_apr(ret, ctl, atm_apr, iqa, ipa, s_a_inv);
  matrix_invert(s_a_inv);

  /* Get measurement errors... */
  set_cov_meas(ret, ctl, obs_meas, sig_noise, sig_formod, sig_eps_inv);

  /* Determine dx = x_i - x_a and dy = y - F(x_i) ... */
  gsl_vector_memcpy(dx, x_i);
  gsl_vector_sub(dx, x_a);
  gsl_vector_memcpy(dy, y_m);
  gsl_vector_sub(dy, y_i);

  /* Compute cost function... */
  *chisq = cost_function(dx, dy, s_a_inv, sig_eps_inv);

  /* Compute initial kernel... */
  kernel(ctl, tbl, atm_i, obs_i, k_i);

  /* ------------------------------------------------------------
     Levenberg-Marquardt minimization...
     ------------------------------------------------------------ */

  /* Outer loop... */
  for (int it = 1; it <= ret->conv_itmax; it++) {

    /* Store current cost function value... */
    double chisq_old = *chisq;

    /* Compute kernel matrix K_i... */
    if (it > 1 && it % ret->kernel_recomp == 0)
      kernel(ctl, tbl, atm_i, obs_i, k_i);

    /* Compute K_i^T * S_eps^{-1} * K_i ... */
    if (it == 1 || it % ret->kernel_recomp == 0)
      matrix_product(k_i, sig_eps_inv, 1, cov);

    /* Determine b = K_i^T * S_eps^{-1} * dy - S_a^{-1} * dx ... */
    for (size_t i = 0; i < m; i++)
      gsl_vector_set(y_aux, i, gsl_vector_get(dy, i)
		     * gsl_pow_2(gsl_vector_get(sig_eps_inv, i)));
    gsl_blas_dgemv(CblasTrans, 1.0, k_i, y_aux, 0.0, b);
    gsl_blas_dgemv(CblasNoTrans, -1.0, s_a_inv, dx, 1.0, b);

    /* Inner loop... */
    for (int it2 = 0; it2 < 20; it2++) {

      /* Compute A = (1 + lmpar) * S_a^{-1} + K_i^T * S_eps^{-1} * K_i ... */
      gsl_matrix_memcpy(a, s_a_inv);
      gsl_matrix_scale(a, 1 + lmpar);
      gsl_matrix_add(a, cov);

      /* Solve A * x_step = b by means of Cholesky decomposition... */
      gsl_linalg_cholesky_decomp(a);
      gsl_linalg_cholesky_solve(a, b, x_step);

      /* Update atmospheric state... */
      gsl_vector_add(x_i, x_step);
      copy_atm(ctl, atm_i, atm_apr, 0);
      copy_obs(ctl, obs_i, obs_meas, 0);
      x2atm(ctl, x_i, atm_i);

      /* Check atmospheric state... */
      for (int ip = 0; ip < atm_i->np; ip++) {
	atm_i->p[ip] = GSL_MIN(GSL_MAX(atm_i->p[ip], 5e-7), 5e4);
	atm_i->t[ip] = GSL_MIN(GSL_MAX(atm_i->t[ip], 100), 400);
	for (int ig = 0; ig < ctl->ng; ig++)
	  atm_i->q[ig][ip] = GSL_MIN(GSL_MAX(atm_i->q[ig][ip], 0), 1);
	for (int iw = 0; iw < ctl->nw; iw++)
	  atm_i->k[iw][ip] = GSL_MAX(atm_i->k[iw][ip], 0);
      }

      /* Forward calculation... */
      formod(ctl, tbl, atm_i, obs_i);
      obs2y(ctl, obs_i, y_i, NULL, NULL);

      /* Determine dx = x_i - x_a and dy = y - F(x_i) ... */
      gsl_vector_memcpy(dx, x_i);
      gsl_vector_sub(dx, x_a);
      gsl_vector_memcpy(dy, y_m);
      gsl_vector_sub(dy, y_i);

      /* Compute cost function... */
      *chisq = cost_function(dx, dy, s_a_inv, sig_eps_inv);

      /* Modify Levenberg-Marquardt parameter... */
      if (*chisq > chisq_old) {
	lmpar *= 10;
	gsl_vector_sub(x_i, x_step);
      } else {
	lmpar /= 10;
	break;
      }
    }

    /* Get normalized step size in state space... */
    gsl_blas_ddot(x_step, b, &disq);
    disq /= (double) n;

    /* Convergence test... */
    if ((it == 1 || it % ret->kernel_recomp == 0) && disq < ret->conv_dmin)
      break;
  }

  /* ------------------------------------------------------------
     Finalize...
     ------------------------------------------------------------ */

  gsl_matrix_free(a);
  gsl_matrix_free(cov);
  gsl_matrix_free(k_i);
  gsl_matrix_free(s_a_inv);

  gsl_vector_free(b);
  gsl_vector_free(dx);
  gsl_vector_free(dy);
  gsl_vector_free(sig_eps_inv);
  gsl_vector_free(sig_formod);
  gsl_vector_free(sig_noise);
  gsl_vector_free(x_a);
  gsl_vector_free(x_i);
  gsl_vector_free(x_step);
  gsl_vector_free(y_aux);
  gsl_vector_free(y_i);
  gsl_vector_free(y_m);
}

/*****************************************************************************/

void read_nc(
  char *filename,
  ncd_t *ncd) {

  int varid;

  /* Open netCDF file... */
  printf("Read netCDF file: %s\n", filename);
  NC(nc_open(filename, NC_WRITE, &ncd->ncid));

  /* Read Level-1 data... */
  NC(nc_inq_varid(ncd->ncid, "l1_time", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l1_time[0]));
  NC(nc_inq_varid(ncd->ncid, "l1_lon", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l1_lon[0]));
  NC(nc_inq_varid(ncd->ncid, "l1_lat", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l1_lat[0]));
  NC(nc_inq_varid(ncd->ncid, "l1_sat_z", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l1_sat_z));
  NC(nc_inq_varid(ncd->ncid, "l1_sat_lon", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l1_sat_lon));
  NC(nc_inq_varid(ncd->ncid, "l1_sat_lat", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l1_sat_lat));
  NC(nc_inq_varid(ncd->ncid, "l1_nu", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l1_nu));
  NC(nc_inq_varid(ncd->ncid, "l1_rad", &varid));
  NC(nc_get_var_float(ncd->ncid, varid, ncd->l1_rad[0][0]));

  /* Read Level-2 data... */
  NC(nc_inq_varid(ncd->ncid, "l2_z", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l2_z[0][0]));
  NC(nc_inq_varid(ncd->ncid, "l2_press", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l2_p));
  NC(nc_inq_varid(ncd->ncid, "l2_temp", &varid));
  NC(nc_get_var_double(ncd->ncid, varid, ncd->l2_t[0][0]));
}

/*****************************************************************************/

void read_ret_ctl(
  int argc,
  char *argv[],
  ctl_t *ctl,
  ret_t *ret) {

  /* Iteration control... */
  ret->kernel_recomp =
    (int) scan_ctl(argc, argv, "KERNEL_RECOMP", -1, "3", NULL);
  ret->conv_itmax = (int) scan_ctl(argc, argv, "CONV_ITMAX", -1, "30", NULL);
  ret->conv_dmin = scan_ctl(argc, argv, "CONV_DMIN", -1, "0.1", NULL);

  for (int id = 0; id < ctl->nd; id++)
    ret->err_formod[id] = scan_ctl(argc, argv, "ERR_FORMOD", id, "0", NULL);

  for (int id = 0; id < ctl->nd; id++)
    ret->err_noise[id] = scan_ctl(argc, argv, "ERR_NOISE", id, "0", NULL);

  ret->err_press = scan_ctl(argc, argv, "ERR_PRESS", -1, "0", NULL);
  ret->err_press_cz = scan_ctl(argc, argv, "ERR_PRESS_CZ", -1, "-999", NULL);
  ret->err_press_ch = scan_ctl(argc, argv, "ERR_PRESS_CH", -1, "-999", NULL);

  ret->err_temp = scan_ctl(argc, argv, "ERR_TEMP", -1, "0", NULL);
  ret->err_temp_cz = scan_ctl(argc, argv, "ERR_TEMP_CZ", -1, "-999", NULL);
  ret->err_temp_ch = scan_ctl(argc, argv, "ERR_TEMP_CH", -1, "-999", NULL);

  for (int ig = 0; ig < ctl->ng; ig++) {
    ret->err_q[ig] = scan_ctl(argc, argv, "ERR_Q", ig, "0", NULL);
    ret->err_q_cz[ig] = scan_ctl(argc, argv, "ERR_Q_CZ", ig, "-999", NULL);
    ret->err_q_ch[ig] = scan_ctl(argc, argv, "ERR_Q_CH", ig, "-999", NULL);
  }

  for (int iw = 0; iw < ctl->nw; iw++) {
    ret->err_k[iw] = scan_ctl(argc, argv, "ERR_K", iw, "0", NULL);
    ret->err_k_cz[iw] = scan_ctl(argc, argv, "ERR_K_CZ", iw, "-999", NULL);
    ret->err_k_ch[iw] = scan_ctl(argc, argv, "ERR_K_CH", iw, "-999", NULL);
  }
}

/*****************************************************************************/

void set_cov_apr(
  ret_t *ret,
  ctl_t *ctl,
  atm_t *atm,
  int *iqa,
  int *ipa,
  gsl_matrix *s_a) {

  /* Get sizes... */
  const size_t n = s_a->size1;

  /* Allocate... */
  gsl_vector *x_a = gsl_vector_alloc(n);

  /* Get sigma vector... */
  atm2x(ctl, atm, x_a, NULL, NULL);
  for (size_t i = 0; i < n; i++) {
    if (iqa[i] == IDXP)
      gsl_vector_set(x_a, i, ret->err_press / 100 * gsl_vector_get(x_a, i));
    if (iqa[i] == IDXT)
      gsl_vector_set(x_a, i, ret->err_temp);
    for (int ig = 0; ig < ctl->ng; ig++)
      if (iqa[i] == IDXQ(ig))
	gsl_vector_set(x_a, i, ret->err_q[ig] / 100 * gsl_vector_get(x_a, i));
    for (int iw = 0; iw < ctl->nw; iw++)
      if (iqa[i] == IDXK(iw))
	gsl_vector_set(x_a, i, ret->err_k[iw]);
  }

  /* Check standard deviations... */
  for (size_t i = 0; i < n; i++)
    if (gsl_pow_2(gsl_vector_get(x_a, i)) <= 0)
      ERRMSG("Check a priori data (zero standard deviation)!");

  /* Initialize diagonal covariance... */
  gsl_matrix_set_zero(s_a);
  for (size_t i = 0; i < n; i++)
    gsl_matrix_set(s_a, i, i, gsl_pow_2(gsl_vector_get(x_a, i)));

  /* Loop over matrix elements... */
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if (i != j && iqa[i] == iqa[j]) {

	/* Initialize... */
	double cz = 0;
	double ch = 0;

	/* Set correlation lengths for pressure... */
	if (iqa[i] == IDXP) {
	  cz = ret->err_press_cz;
	  ch = ret->err_press_ch;
	}

	/* Set correlation lengths for temperature... */
	if (iqa[i] == IDXT) {
	  cz = ret->err_temp_cz;
	  ch = ret->err_temp_ch;
	}

	/* Set correlation lengths for volume mixing ratios... */
	for (int ig = 0; ig < ctl->ng; ig++)
	  if (iqa[i] == IDXQ(ig)) {
	    cz = ret->err_q_cz[ig];
	    ch = ret->err_q_ch[ig];
	  }

	/* Set correlation lengths for extinction... */
	for (int iw = 0; iw < ctl->nw; iw++)
	  if (iqa[i] == IDXK(iw)) {
	    cz = ret->err_k_cz[iw];
	    ch = ret->err_k_ch[iw];
	  }

	/* Compute correlations... */
	if (cz > 0 && ch > 0) {

	  /* Get Cartesian coordinates... */
	  double x0[3], x1[3];
	  geo2cart(0, atm->lon[ipa[i]], atm->lat[ipa[i]], x0);
	  geo2cart(0, atm->lon[ipa[j]], atm->lat[ipa[j]], x1);

	  /* Compute correlations... */
	  const double rho =
	    exp(-DIST(x0, x1) / ch -
		fabs(atm->z[ipa[i]] - atm->z[ipa[j]]) / cz);

	  /* Set covariance... */
	  gsl_matrix_set(s_a, i, j, gsl_vector_get(x_a, i)
			 * gsl_vector_get(x_a, j) * rho);
	}
      }

  /* Free... */
  gsl_vector_free(x_a);
}

/*****************************************************************************/

void set_cov_meas(
  ret_t *ret,
  ctl_t *ctl,
  obs_t *obs,
  gsl_vector *sig_noise,
  gsl_vector *sig_formod,
  gsl_vector *sig_eps_inv) {

  static obs_t obs_err;

  /* Get size... */
  const size_t m = sig_eps_inv->size;

  /* Noise error (always considered in retrieval fit)... */
  copy_obs(ctl, &obs_err, obs, 1);
  for (int ir = 0; ir < obs_err.nr; ir++)
    for (int id = 0; id < ctl->nd; id++)
      obs_err.rad[id][ir]
	= (gsl_finite(obs->rad[id][ir]) ? ret->err_noise[id] : GSL_NAN);
  obs2y(ctl, &obs_err, sig_noise, NULL, NULL);

  /* Forward model error (always considered in retrieval fit)... */
  copy_obs(ctl, &obs_err, obs, 1);
  for (int ir = 0; ir < obs_err.nr; ir++)
    for (int id = 0; id < ctl->nd; id++)
      obs_err.rad[id][ir]
	= fabs(ret->err_formod[id] / 100 * obs->rad[id][ir]);
  obs2y(ctl, &obs_err, sig_formod, NULL, NULL);

  /* Total error... */
  for (size_t i = 0; i < m; i++)
    gsl_vector_set(sig_eps_inv, i,
		   1 / sqrt(gsl_pow_2(gsl_vector_get(sig_noise, i))
			    + gsl_pow_2(gsl_vector_get(sig_formod, i))));

  /* Check standard deviations... */
  for (size_t i = 0; i < m; i++)
    if (gsl_vector_get(sig_eps_inv, i) <= 0)
      ERRMSG("Check measurement errors (zero standard deviation)!");
}

/*****************************************************************************/

void write_nc(
  char *filename,
  ncd_t *ncd) {

  int dimid[10], p_id, t_id, z_id;

  /* Create netCDF file... */
  printf("Write netCDF file: %s\n", filename);

  /* Read existing dimensions... */
  NC(nc_inq_dimid(ncd->ncid, "L1_NTRACK", &dimid[0]));
  NC(nc_inq_dimid(ncd->ncid, "L1_NXTRACK", &dimid[1]));

  /* Set define mode... */
  NC(nc_redef(ncd->ncid));

  /* Set new dimensions... */
  if (nc_inq_dimid(ncd->ncid, "RET_NP", &dimid[2]) != NC_NOERR)
    NC(nc_def_dim(ncd->ncid, "RET_NP", (size_t) ncd->np, &dimid[2]));

  /* Set new variables... */
  add_var(ncd->ncid, "ret_z", "km", "altitude", NC_FLOAT, &dimid[2], &z_id,
	  1);
  add_var(ncd->ncid, "ret_press", "hPa", "pressure", NC_FLOAT, dimid, &p_id,
	  2);
  add_var(ncd->ncid, "ret_temp", "K", "temperature", NC_FLOAT, dimid, &t_id,
	  3);

  /* Leave define mode... */
  NC(nc_enddef(ncd->ncid));

  /* Write data... */
  NC(nc_put_var_float(ncd->ncid, z_id, ncd->ret_z));
  NC(nc_put_var_float(ncd->ncid, p_id, ncd->ret_p));
  NC(nc_put_var_float(ncd->ncid, t_id, ncd->ret_t));

  /* Close netCDF file... */
  NC(nc_close(ncd->ncid));
}
