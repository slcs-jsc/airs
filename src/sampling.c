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
  Estimate sampling patterns of AIRS.
*/

#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  double d, dmin, dmax, dmu, x0[3], x1[3], x2[3];

  int i, itrack, ixtrack, n;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <pert.nc>");

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], "4mu", pert);

  /* Init... */
  dmin = 1e100;
  dmax = -1e100;
  dmu = 0;
  n = 0;

  /* Get swath width... */
  for (itrack = 0; itrack < pert->ntrack; itrack++) {
    geo2cart(0, pert->lon[itrack][0], pert->lat[itrack][0], x0);
    geo2cart(0, pert->lon[itrack][pert->nxtrack - 1],
	     pert->lat[itrack][pert->nxtrack - 1], x1);
    d = 2. * RE * asin(DIST(x0, x1) / (2. * RE));
    dmin = GSL_MIN(dmin, d);
    dmax = GSL_MAX(dmax, d);
    dmu += d;
    n++;
  }

  /* Write output... */
  printf("\nmean_swath_width=    %.1f km\n", dmu / n);
  printf("minimum_swath_width= %.1f km\n", dmin);
  printf("maximum_swath_width= %.1f km\n", dmax);

  /* Init... */
  dmin = 1e100;
  dmax = -1e100;
  dmu = 0;
  n = 0;

  /* Get across-track sampling distances... */
  for (itrack = 0; itrack < pert->ntrack; itrack++) {
    for (ixtrack = 0; ixtrack < pert->nxtrack - 1; ixtrack++) {
      geo2cart(0, pert->lon[itrack][ixtrack], pert->lat[itrack][ixtrack], x0);
      geo2cart(0, pert->lon[itrack][ixtrack + 1],
	       pert->lat[itrack][ixtrack + 1], x1);
      d = 2. * RE * asin(DIST(x0, x1) / (2. * RE));
      dmin = GSL_MIN(dmin, d);
      dmax = GSL_MAX(dmax, d);
      dmu += d;
      n++;
    }
  }

  /* Write output... */
  printf("\nmean_across_track_sampling_distance=    %.1f km\n", dmu / n);
  printf("minimum_across_track_sampling_distance= %.1f km\n", dmin);
  printf("maximum_across_track_sampling_distance= %.1f km\n", dmax);

  /* Init... */
  dmin = 1e100;
  dmax = -1e100;
  dmu = 0;
  n = 0;

  /* Get along-track sampling distances... */
  for (itrack = 0; itrack < pert->ntrack - 1; itrack++) {
    for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {
      geo2cart(0, pert->lon[itrack][ixtrack], pert->lat[itrack][ixtrack], x0);
      geo2cart(0, pert->lon[itrack + 1][ixtrack],
	       pert->lat[itrack + 1][ixtrack], x1);
      d = 2. * RE * asin(DIST(x0, x1) / (2. * RE));
      dmin = GSL_MIN(dmin, d);
      dmax = GSL_MAX(dmax, d);
      dmu += d;
      n++;
    }
  }

  /* Write output... */
  printf("\nmean_along_track_sampling_distance=    %.1f km\n", dmu / n);
  printf("minimum_along_track_sampling_distance= %.1f km\n", dmin);
  printf("maximum_along_track_sampling_distance= %.1f km\n", dmax);

  /* Init... */
  dmin = 1e100;
  dmax = -1e100;
  dmu = 0;
  n = 0;

  /* Get angle between along-track and across-track direction... */
  for (itrack = 0; itrack < pert->ntrack - 1; itrack++) {
    geo2cart(0, pert->lon[itrack][pert->nxtrack / 2],
	     pert->lat[itrack][pert->nxtrack / 2], x0);
    geo2cart(0, pert->lon[itrack][pert->nxtrack / 2 + 1],
	     pert->lat[itrack][pert->nxtrack / 2 + 1], x1);
    geo2cart(0, pert->lon[itrack + 1][pert->nxtrack / 2],
	     pert->lat[itrack + 1][pert->nxtrack / 2], x2);
    for (i = 0; i < 3; i++) {
      x1[i] -= x0[i];
      x2[i] -= x0[i];
    }
    d = acos(DOTP(x1, x2) / (NORM(x1) * NORM(x2))) * 180. / M_PI;
    dmin = GSL_MIN(dmin, d);
    dmax = GSL_MAX(dmax, d);
    dmu += d;
    n++;
  }

  /* Write output... */
  printf("\nmean_across_track_angle=    %.1f deg\n", dmu / n);
  printf("minimum_across_track_angle= %.1f deg\n", dmin);
  printf("maximum_across_track_angle= %.1f deg\n", dmax);

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}
