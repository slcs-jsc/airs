# AGENTS.md

This file provides repository-specific guidance for coding agents working on the AIRS Code Collection.

## Project Identity

- Project name: `AIRS Code Collection`
- Language: C
- Build system: `make`
- Main source tree: `src/`
- Documentation trees:
  - MkDocs source: `docs/manual/`
  - MkDocs config: `docs/mkdocs.yml`
  - Doxygen config: `docs/Doxyfile`

## Repository Layout

- `src/`: executables, `libairs.*`, and the copied `jurassic.*` dependency sources
- `libs/`: bundled third-party archives plus `build.sh` to build local dependencies
- `tests/`: regression tests and sample AIRS input data
- `docs/`: MkDocs and Doxygen configuration plus generated documentation output

## Build And Test

Run commands from the repository root unless stated otherwise.

- Build bundled libraries if needed:
```bash
cd libs
./build.sh
```

- Build the code:
```bash
cd src
make
```

- Run regression tests:
```bash
cd src
make check
```

- Build MkDocs manual:
```bash
cd src
make mkdocs
```

- Build Doxygen manual:
```bash
cd src
make doxygen
```

## Important Build Notes

- The default `Makefile` uses strict warnings and `-Werror`. Small warning regressions will fail the build.
- Default include and library search paths assume locally built libraries under `libs/build/` plus standard Debian or Ubuntu library locations.
- `nlte` and `retrieval` are built with `mpicc.openmpi`.
- Test scripts expect:
```bash
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
```

## Documentation Conventions

- Keep the project name as `AIRS Code Collection`.
- Use the existing MkDocs content under `docs/manual/`; do not create a separate docs tree unless explicitly requested.
- The Doxygen mainpage is defined in `src/libairs.h`.
- `jurassic.h` and `jurassic.c` are copied into this repository for implementation purposes, but Doxygen should treat JURASSIC as an external upstream dependency.

## JURASSIC Dependency

- `src/jurassic.h` and `src/jurassic.c` come from the upstream JURASSIC project.
- Prefer not to make incidental edits to the copied JURASSIC sources unless the task specifically requires it.
- When describing the dependency in documentation, refer to JURASSIC as an upstream radiative transfer model and retrieval framework.

## Editing Guidance

- Prefer targeted edits over broad refactors.
- Preserve the existing C style and file header structure.
- Do not silently rewrite large sets of source files unless the task explicitly calls for normalization.
- If changing generated docs behavior, verify with the relevant build command rather than assuming the config is correct.

## Verification Expectations

For documentation changes:

- run `make mkdocs` for MkDocs changes
- run `make doxygen` for Doxygen changes

For source changes:

- run at least the most relevant target from `make` or `make check` when feasible

## Generated Output

- MkDocs output is generated in `docs/site/`
- Doxygen output is generated in `docs/docs/html/`

Do not hand-edit generated HTML unless the task explicitly asks for generated artifacts to be modified directly.
