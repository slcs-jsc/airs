# Command Reference

This page groups the executables listed in `src/Makefile`. It is a user-level summary rather than a full API description.

## AIRS Radiance And Footprint Tools

- `spec2tab`: extract radiance spectra from a Level-1B granule. Usage: `<airs_l1b_file> [index <track> <xtrack> | geo <lon> <lat>] <spec.tab>`.
- `spec_qual`: extract quality flags for one AIRS footprint. Usage: `<airs_l1b_file> <track> <xtrack> <qual.tab>`.
- `map_rad`: extract radiance maps from one or two Level-1B granules. Usage: `<ctl> <l1b_file1> <l1b_file2> <nu> <wave.tab>`.
- `orbit`: extract AIRS/Aqua orbit information. Usage: `<orbit.tab> <airs_l1b_file> [<airs_l1b_file2> ...]`.
- `optimize_btd`: optimize brightness-temperature differences.
- `volcano`: detect volcanic ash and sulfur dioxide from Level-1B radiances. Usage: `<out.tab> <l1b_file1> [<l1b_file2> ...]`.

## Perturbation, Noise, And Wave Analysis

- `perturbation`: create perturbation NetCDF products from Level-1B granules. Usage: `<out.nc> <l1b_file1> [<l1b_file2> ...]`.
- `map_pert`: extract perturbation maps. Usage: `<ctl> <pert.nc> <map.tab>`.
- `noise_pert`: estimate noise from perturbation products. Usage: `<ctl> <pert.nc> <noise.tab>`.
- `variance`: calculate brightness-temperature variances. Usage: `<ctl> <var.tab> <pert1.nc> [<pert2.nc> ...]`.
- `events`: identify gravity-wave events. Usage: `<ctl> <events.tab> <pert1.nc> [<pert2.nc> ...]`.
- `spec_ana`: run spectral analysis of gravity-wave perturbations. Usage: `<ctl> <pert.nc> <spec.tab>`.
- `spec_synth`: run the same spectral analysis workflow on synthetic data. Usage: `<ctl> <spec.tab>`.
- `sampling`: estimate AIRS sampling patterns. Usage: `<ctl> <pert.nc>`.
- `var1d`: estimate horizontal wavelength sensitivity. Usage: `<width> <n> <lxmin> <lxmax> <dlx> <fwhm> <dim>`.
- `rayt`: run the 2-D gravity-wave ray-tracing tool.

## Retrieval Preparation And Post-Processing

- `extract`: extract radiance data for retrieval inputs. Usage: `<airs_l1_file> <airs_l2_file> <out.nc>`.
- `retrieval`: AIRS retrieval processor. Usage: `<ctl> <filelist>`.
- `nlte`: retrieve non-LTE index. Usage: `<ctl> <filelist>`.
- `ret2tab`: write retrieval data to ASCII. Usage: `<airs_l2_file> <layer> <airs.tab>`.
- `map_ret`: extract retrieval maps. Usage: `<ctl> <airs.nc> <map.tab>`.
- `noise_ret`: estimate noise from retrieval data. Usage: `<ctl> <airs.nc> <noise.tab>`.
- `diff_apr`: compare retrieval and a priori data. Usage: `<ctl> <airs.nc> <airs2.nc> <diff.tab>`.
- `diff_ret`: compare retrieval outputs. Usage: `<ctl> <airs.nc> <airs2.nc> <diff.tab>`.
- `select_ret`: find retrieval results for a geographic selection. Usage: `<ctl> [<airs1.nc> <airs2.nc> ...]`.
- `zm_ret`: calculate zonal means for retrieval data. Usage: `<ctl> <zm.tab> <airs1.nc> [<airs2.nc> ...]`.

## Case-Specific Analysis Programs

- `hurricane`: analyze gravity-wave data for tropical cyclones.
- `island`: analyze gravity-wave data for remote islands.
- `overpass`: find AIRS/Aqua overpasses for a target location. Usage: `<ctl> <pert.nc> <lon0> <lat0> <overpass.tab>`.
- `issifm`: simulate AIRS observations from atmospheric model output.

## Utility Programs

- `day2doy`: convert calendar date to day-of-year.
- `doy2day`: convert day-of-year to calendar date.
- `time2jsec`: convert date and time to Julian seconds.
- `jsec2time`: convert Julian seconds to date and time.
- `distance`: calculate distance between two geolocations.
- `sza`: calculate solar zenith angle.

## Shared Libraries

- `libairs.c` and `libairs.h`: AIRS-specific data structures and helper routines.
- `jurassic.c` and `jurassic.h`: retrieval framework used by `retrieval` and `nlte`.

For implementation-level details, use the Doxygen output linked on the [Links](links.md) page.
