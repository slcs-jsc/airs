# Testing

## Available Regression Tests

The source `Makefile` defines two regression targets:

- `pert_test`
- `rad_test`

Run both with:

```bash
cd src
make check
```

Or run one test directly:

```bash
cd tests/rad_test
./run.sh
```

## Test Data

Sample AIRS Level-1B granules used by the tests are stored in `tests/data`.

Reference outputs are stored under:

- `tests/rad_test/data.ref`
- `tests/pert_test/data.ref`

## Test Environment

The regression scripts set:

```bash
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
```

That is the expected runtime pattern when relying on locally built libraries from `libs/build`.

## What The Tests Cover

- `rad_test`: spectrum extraction, quality extraction, radiance mapping, and orbit extraction.
- `pert_test`: perturbation generation, perturbation mapping, noise estimation, variance analysis, and event detection.

## Documentation CI

The repository workflow `.github/workflows/docs.yml` installs `mkdocs` and `mkdocs-material`, builds the MkDocs site, builds the Doxygen output, and publishes both together via GitHub Pages.
