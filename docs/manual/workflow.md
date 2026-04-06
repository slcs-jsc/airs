# Workflow

This page summarizes the practical workflows already encoded in the repository test suite.

## Radiance Inspection Workflow

The `rad_test` regression demonstrates a minimal Level-1B radiance workflow:

1. Extract a spectrum with `spec2tab`.
2. Extract quality flags with `spec_qual`.
3. Build a radiance map with `map_rad`.
4. Extract orbit information with `orbit`.

Representative commands:

```bash
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

../../src/spec2tab ../data/AIRS...02637.hdf index 60 44 data/spec.tab
../../src/spec_qual ../data/AIRS...02637.hdf 60 44 data/qual.tab
../../src/map_rad - ../data/AIRS...02637.hdf ../data/AIRS...02639.hdf 2338.43 data/wave.tab
../../src/orbit data/orbit.tab ../data/AIRS...02637.hdf
```

The sample granules live in `tests/data`.

## Perturbation And Gravity-Wave Workflow

The `pert_test` regression demonstrates the gravity-wave analysis chain:

1. Create a perturbation product with `perturbation`.
2. Extract perturbation maps with `map_pert`.
3. Estimate noise with `noise_pert`.
4. Estimate variance with `variance`.
5. Detect events with `events`.

Representative commands:

```bash
../../src/perturbation data/pert.nc ../data/AIRS...02637.hdf ../data/AIRS...02639.hdf
../../src/map_pert - data/pert.nc data/map_4mu.tab PERTNAME 4mu
../../src/noise_pert - data/pert.nc data/noise_4mu.tab PERTNAME 4mu
../../src/variance - data/var_4mu.tab data/pert.nc PERTNAME 4mu NX 60 NY 30
../../src/events - data/events_4mu.tab data/pert.nc PERTNAME 4mu VARMIN 0.2
```

The same workflow is used for these predefined perturbation channel sets:

- `4mu`
- `15mu_low`
- `15mu_high`

## Retrieval-Oriented Workflow

The retrieval side of the code base is centered on:

- `extract` for preparing AIRS data for retrievals
- `retrieval` for the main AIRS retrieval processor
- `nlte` for non-LTE index retrieval
- post-processing tools such as `ret2tab`, `map_ret`, `noise_ret`, `diff_apr`, `diff_ret`, `select_ret`, and `zm_ret`

In practice, start by preparing a control file and a file list, then run `retrieval` or `nlte` before using the downstream diagnostics.

## Runtime Notes

- Most analysis commands accept `-` as the first argument instead of a control file, allowing parameters to be supplied directly on the command line.
- Several programs use OpenMP; setting `OMP_NUM_THREADS` is recommended for reproducible runs.
- The test scripts prepend `../../libs/build/lib` to `LD_LIBRARY_PATH`, which is the expected pattern when relying on locally built libraries.
