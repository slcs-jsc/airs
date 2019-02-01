#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static ret_t ret, ret2;

  static FILE *out;

  static double mean[NPG], sigma[NPG], min[NPG], max[NPG],
    tt[NPG], lon[NPG], lat[NPG], temp[NPG], press[NPG];

  static int ids, ip;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <airs.nc> <airs2.nc> <diff.tab>");

  /* Read AIRS data... */
  read_retr(argv[2], &ret);
  read_retr(argv[3], &ret2);

  /* Compute differences... */
  for (ids = 0; ids < ret.nds; ids++)
    for (ip = 0; ip < ret.np; ip++) {
      if (ret.time[ids][ip] != ret2.time[ids][ip] ||
	  ret.lon[ids][ip] != ret2.lon[ids][ip] ||
	  ret.lat[ids][ip] != ret2.lat[ids][ip])
	ERRMSG("Data files do not match!");
      tt[ip] += ret.time[ids][ip];
      lon[ip] += ret.lon[ids][ip];
      lat[ip] += ret.lat[ids][ip];
      press[ip] += ret.p[ids][ip];
      temp[ip] += ret.t[ids][ip];
      mean[ip] += ret2.t[ids][ip] - ret.t[ids][ip];
      sigma[ip] += gsl_pow_2(ret2.t[ids][ip] - ret.t[ids][ip]);
      min[ip] = GSL_MIN(min[ip], ret2.t[ids][ip] - ret.t[ids][ip]);
      max[ip] = GSL_MAX(max[ip], ret2.t[ids][ip] - ret.t[ids][ip]);
    }

  /* Create output file... */
  printf("Write retrieval differences data: %s\n", argv[4]);
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
  for (ip = 0; ip < ret.np; ip++)
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g\n",
	    tt[ip] / ret.nds, ret.z[0][ip], lon[ip] / ret.nds,
	    lat[ip] / ret.nds, press[ip] / ret.nds, temp[ip] / ret.nds,
	    mean[ip] / ret.nds,
	    sqrt(sigma[ip] / ret.nds - gsl_pow_2(mean[ip] / ret.nds)),
	    min[ip], max[ip]);

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
