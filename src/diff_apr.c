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
  Extract differences between retrieval and a priori data.
*/

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

/* Number of AIRS radiance channels (don't change). */
#define L1_NCHAN 34

/* Along-track size of AIRS radiance granule (don't change). */
#define L1_NTRACK 135

/* Across-track size of AIRS radiance granule (don't change). */
#define L1_NXTRACK 90

/* Number of AIRS pressure layers (don't change). */
#define L2_NLAY 27

/* Along-track size of AIRS retrieval granule (don't change). */
#define L2_NTRACK 45

/* Across-track size of AIRS retrieval granule (don't change). */
#define L2_NXTRACK 30

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/* Buffer for netCDF data. */
typedef struct {

  /* NetCDF file ID. */
  int ncid;

  /* Number of retrieval altitudes. */
  int np;

  /* Time (seconds since 2000-01-01T00:00Z). */
  double l1_time[L1_NTRACK][L1_NXTRACK];

  /* Footprint longitude [deg]. */
  double l1_lon[L1_NTRACK][L1_NXTRACK];

  /* Footprint latitude [deg]. */
  double l1_lat[L1_NTRACK][L1_NXTRACK];

  /* Satellite altitude [km]. */
  double l1_sat_z[L1_NTRACK];

  /* Satellite longitude [deg]. */
  double l1_sat_lon[L1_NTRACK];

  /* Satellite latitude [deg]. */
  double l1_sat_lat[L1_NTRACK];

  /* Channel frequencies [cm^-1]. */
  double l1_nu[L1_NCHAN];

  /* Radiance [W/(m^2 sr cm^-1)]. */
  float l1_rad[L1_NTRACK][L1_NXTRACK][L1_NCHAN];

  /* Altitude [km]. */
  double l2_z[L2_NTRACK][L2_NXTRACK][L2_NLAY];

  /* Pressure [hPa]. */
  double l2_p[L2_NLAY];

  /* Temperature [K]. */
  double l2_t[L2_NTRACK][L2_NXTRACK][L2_NLAY];

  /* Altitude [km]. */
  float ret_z[NP];

  /* Pressure [hPa]. */
  float ret_p[L1_NTRACK * L1_NXTRACK];

  /* Temperature [K]. */
  float ret_t[L1_NTRACK * L1_NXTRACK * NP];

} ncd_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Read netCDF file. */
void read_nc(
  char *filename,
  ncd_t * ncd);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;

  static ncd_t ncd, ncd2;

  static FILE *out;

  static double mean[L2_NLAY], sigma[L2_NLAY], min[L2_NLAY], max[L2_NLAY],
    tt[L2_NLAY], lon[L2_NLAY], lat[L2_NLAY], temp[L2_NLAY], press[L2_NLAY],
    z[L2_NLAY], tip;

  static int idx, ip, itrack, ixtrack;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <airs.nc> <airs2.nc> <diff.tab>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read netCDF files... */
  read_nc(argv[2], &ncd);
  read_nc(argv[3], &ncd2);

  /* Compute differences... */
  for (itrack = 0; itrack < L2_NTRACK; itrack++)
    for (ixtrack = 0; ixtrack < L2_NXTRACK; ixtrack++) {
      for (ip = 0; ip < L2_NLAY; ip++) {
	if (ncd.l1_time[3 * itrack + 1][3 * ixtrack + 1] !=
	    ncd2.l1_time[3 * itrack + 1][3 * ixtrack + 1]
	    || ncd.l1_lon[3 * itrack + 1][3 * ixtrack + 1] !=
	    ncd2.l1_lon[3 * itrack + 1][3 * ixtrack + 1]
	    || ncd.l1_lat[3 * itrack + 1][3 * ixtrack + 1] !=
	    ncd2.l1_lat[3 * itrack + 1][3 * ixtrack + 1])
	  ERRMSG("Data files do not match!");
	tt[ip] += ncd.l1_time[3 * itrack + 1][3 * ixtrack + 1];
	lon[ip] += ncd.l1_lon[3 * itrack + 1][3 * ixtrack + 1];
	lat[ip] += ncd.l1_lat[3 * itrack + 1][3 * ixtrack + 1];
	z[ip] += ncd.l2_z[itrack][ixtrack][ip];
	press[ip] += ncd.l2_p[ip];
	temp[ip] += ncd.l2_t[itrack][ixtrack][ip];
	idx =
	  locate_irr(ncd2.l2_z[itrack][ixtrack], L2_NLAY,
		     ncd.l2_z[itrack][ixtrack][ip]);
	tip =
	  LIN(ncd2.l2_z[itrack][ixtrack][idx],
	      ncd2.l2_t[itrack][ixtrack][idx],
	      ncd2.l2_z[itrack][ixtrack][idx + 1],
	      ncd2.l2_t[itrack][ixtrack][idx + 1],
	      ncd.l2_z[itrack][ixtrack][ip]);
	mean[ip] += tip - ncd.l2_t[itrack][ixtrack][ip];
	sigma[ip] += gsl_pow_2(tip - ncd.l2_t[itrack][ixtrack][ip]);
	min[ip] = GSL_MIN(min[ip], tip - ncd.l2_t[itrack][ixtrack][ip]);
	max[ip] = GSL_MAX(max[ip], tip - ncd.l2_t[itrack][ixtrack][ip]);
      }
    }

  /* Create output file... */
  printf("Write a priori differences data: %s\n", argv[4]);
  if (!(out = fopen(argv[4], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = pressure (set 1) [hPa]\n"
	  "# $6 = temperature (set 1) [K]\n"
	  "# $7 = temperature difference (mean, set 2 - set 1) [K]\n"
	  "# $8 = temperature difference (sigma, set 2 - set 1) [K]\n"
	  "# $9 = temperature difference (minimum, set 2 - set 1) [K]\n"
	  "# $10 = temperature difference (maximum, set 2 - set 1) [K]\n\n");

  /* Write output... */
  for (ip = 0; ip < L2_NLAY; ip++)
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g\n",
	    tt[ip] / (L2_NTRACK * L2_NXTRACK),
	    z[ip] / (L2_NTRACK * L2_NXTRACK),
	    lon[ip] / (L2_NTRACK * L2_NXTRACK),
	    lat[ip] / (L2_NTRACK * L2_NXTRACK),
	    press[ip] / (L2_NTRACK * L2_NXTRACK),
	    temp[ip] / (L2_NTRACK * L2_NXTRACK),
	    mean[ip] / (L2_NTRACK * L2_NXTRACK),
	    sqrt(sigma[ip] / (L2_NTRACK * L2_NXTRACK) -
		 gsl_pow_2(mean[ip] / (L2_NTRACK * L2_NXTRACK))), min[ip],
	    max[ip]);

  /* Close file... */
  fclose(out);
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
