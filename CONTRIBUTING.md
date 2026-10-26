# Contributing to FibrTraQ

FibrTraQ is a scientific R package for voxel-resolved focal-ROI target profiling of diffusion tractograms. Contributions should preserve the package's current scientific semantics unless a change is explicitly proposed, tested, and documented.

## License

Contributions to FibrTraQ are distributed under the
[GNU General Public License, version 2 or later](LICENSE.md).
Preserve third-party copyright and license notices, and identify any
third-party material included in a contribution.

## Development Setup

1. Install R 4.0 or newer.
2. Install package dependencies from `DESCRIPTION`.
3. From the repository root, load the package with:

```r
devtools::load_all(".")
```

## Before Submitting Changes

Run the focused checks:

```r
devtools::document()
devtools::test()
Rscript inst/scripts/fibrtraq-path-audit.R
```

For release-facing changes, also run:

```r
R CMD build .
R CMD check --as-cran FibrTraQ_*.tar.gz
```

## Scientific Semantics

The primary analysis mode is `all_focal_voxels` with `sample_points` voxelization and full-streamline target-ROI sampled-point intersections. The legacy `first_intersection` mode is retained for sensitivity analyses. Any proposed change to counting, coordinate handling, local/fate semantics, or target-intersection behavior must include regression tests and migration notes.

## Pull Requests

Please include:

- a concise description of the change;
- the tests or validation commands you ran;
- any expected changes to outputs;
- whether the change affects public API, file formats, or manuscript validation artifacts.
