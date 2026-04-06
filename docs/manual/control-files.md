# Control Files

Many AIRS executables use the shared `scan_ctl(...)` parser implemented in `src/jurassic.c`. The first positional argument is typically a control file path named here as `<ctl>`.

## Two Supported Modes

### Use A Control File

If the first argument is not `-`, the program opens that file and reads settings from it.

Example:

```bash
map_pert my_case.ctl data/pert.nc data/map.tab
```

### Use Inline Key/Value Arguments

If the first argument is `-`, no control file is opened and settings are read from the command line instead.

Example:

```bash
map_pert - data/pert.nc data/map.tab PERTNAME 4mu BG_POLY_X 0 BG_POLY_Y 0
```

## Expected File Format

The parser looks for lines with three whitespace-separated fields:

```text
VARIABLE = VALUE
```

The middle field is parsed but otherwise ignored, so the usual convention is `name = value`.

Examples:

```text
PERTNAME = 4mu
BG_POLY_X = 0
BG_POLY_Y = 0
VAR_DH = 100.0
```

## Array Syntax

Some retrieval settings use indexed variables. The parser accepts both:

- `NAME[3] = value`
- `NAME[*] = value`

The wildcard form applies the same value to all indices that query that setting.

## Override Behavior

Command-line values override values from the control file when both are present.

This is useful for parameter sweeps:

```bash
variance base.ctl out.tab pert.nc NX 60 NY 30
```

## Missing Variables

- If a program provides a default value in code, that default is used.
- If no default exists, execution stops with a missing-variable error.

## Common Keys

The exact keys depend on the executable, but frequently used options include:

- `PERTNAME`
- `SET`
- `BG_POLY_X`
- `BG_POLY_Y`
- `BG_SMOOTH_X`
- `BG_SMOOTH_Y`
- `GAUSS_FWHM`
- `VAR_DH`
- `NU`
- `TRACK_MIN`, `TRACK_MAX`
- `XTRACK_MIN`, `XTRACK_MAX`

The best way to discover program-specific keys is to inspect the `scan_ctl(...)` calls near the top of each program in `src/`.
