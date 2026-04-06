# AIRS Code Collection

The AIRS Code Collection provides command-line tools and supporting libraries for processing, analyzing, and retrieving information from NASA Atmospheric InfraRed Sounder observations.

This manual is the user-facing entry point for the repository. It complements the generated Doxygen reference with practical guidance on building the code, understanding the executable layout, and running the main AIRS workflows.

## Repository Overview

- `src/`: C source files, shared libraries, and executable targets.
- `libs/`: bundled dependency archives and the helper script to build local libraries.
- `tests/`: regression workflows with sample AIRS data and reference outputs.
- `docs/`: MkDocs and Doxygen configuration.

## Main Capabilities

- Read AIRS Level-1B radiance granules and derive spectra, quality flags, maps, and orbit information.
- Create perturbation products for gravity-wave analysis.
- Estimate noise, variances, wave events, and spectral characteristics.
- Run AIRS retrieval and non-LTE processing through the shared JURASSIC infrastructure.
- Provide small utility programs for time conversion, geometry, and diagnostics.

## Documentation Map

- [Getting Started](getting-started.md): prerequisites, dependency build, and compilation.
- [Workflow](workflow.md): practical end-to-end usage patterns derived from the repository tests.
- [Control Files](control-files.md): how AIRS tools read `ctl` settings.
- [Command Reference](command-reference.md): grouped summary of the executables built from `src/`.
- [Testing](testing.md): available regression checks and expected inputs.
- [Links](links.md): repository, citation, Doxygen, and related references.

## Contact

Dr. Lars Hoffmann, <l.hoffmann@fz-juelich.de>

Jülich Supercomputing Centre, Forschungszentrum Jülich
