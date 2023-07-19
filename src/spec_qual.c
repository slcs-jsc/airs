#include "libairs.h"

int main(
  int argc,
  char *argv[]) {

  static airs_rad_gran_t airs_rad_gran;

  FILE *out;

  /* Check arguments... */
  if (argc != 5)
    ERRMSG("Give parameters: <airs_l1b_file> <track> <xtrack> <qual.tab>");

  /* Read AIRS data... */
  printf("Read AIRS Level-1B data file: %s\n", argv[1]);
  airs_rad_rdr(argv[1], &airs_rad_gran);

  /* Get indices... */
  int track = atoi(argv[2]);
  int xtrack = atoi(argv[3]);

  /* Check indices... */
  if (track < 0 || track >= AIRS_RAD_GEOTRACK)
    ERRMSG("Along-track index out of range!");
  if (xtrack < 0 || xtrack >= AIRS_RAD_GEOXTRACK)
    ERRMSG("Across-track index out of range!");

  /* Create file... */
  printf("Write quality data file: %s\n", argv[4]);
  if (!(out = fopen(argv[4], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = satellite longitude [deg]\n"
	  "# $3 = satellite latitude [deg]\n"
	  "# $4 = footprint longitude [deg]\n"
	  "# $5 = footprint latitude [deg]\n"
	  "# $6 = wavenumber [cm^-1]\n"
	  "# $7 = radiance [W/(m^2 sr cm^-1)]\n"
	  "# $8 = channel number\n"
	  "# $9 = state\n"
	  "# $10 = excluded channels\n"
	  "# $11 = calibration channel summary\n"
	  "# $12 = calibration flag\n"
	  "# $13 = combined quality flag\n\n");
  
  /* Write data... */
  for (int ichan = 0; ichan < AIRS_RAD_CHANNEL; ichan++) {
    int flag = (airs_rad_gran.state[track][xtrack] != 0)
      || (airs_rad_gran.ExcludedChans[ichan] > 2)
      || (airs_rad_gran.CalChanSummary[ichan] & 8)
      || (airs_rad_gran.CalChanSummary[ichan] & (32 + 64))
      || (airs_rad_gran.CalFlag[track][ichan] & 16);
    fprintf(out, "%.2f %g %g %g %g %g %g %d %d %d %d %d %d\n",
	    airs_rad_gran.Time[track][xtrack] - 220838400,
	    airs_rad_gran.sat_lon[track],
	    airs_rad_gran.sat_lat[track],
	    airs_rad_gran.Longitude[track][xtrack],
	    airs_rad_gran.Latitude[track][xtrack],
	    airs_rad_gran.nominal_freq[ichan],
	    airs_rad_gran.radiances[track][xtrack][ichan] * 1e-3,
	    ichan,
	    airs_rad_gran.state[track][xtrack],
	    airs_rad_gran.ExcludedChans[ichan],
	    airs_rad_gran.CalChanSummary[ichan],
	    airs_rad_gran.CalFlag[track][ichan],
	    flag);
  }
  
  /* Close file... */
  fclose(out);
  
  return EXIT_SUCCESS;
}
