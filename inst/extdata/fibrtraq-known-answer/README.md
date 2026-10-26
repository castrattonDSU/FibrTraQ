# FibrTraQ Known-Answer Example

This package-authored synthetic dataset provides a small, immutable
end-to-end example for FibrTraQ. It contains 14 TrackVis streamlines on
a 9 x 9 x 9 identity-affine voxel grid, a two-voxel focal mask, cortex
and brainstem target masks, a local mask, an endpoint context mask, and
a reference image.

The dataset is designed to exercise focal-voxel assignment, full-path
target intersections, mutually exclusive fate categories, exact and
one-voxel endpoint classes, target penetration, and report generation.
It is not anatomical MRI and must not be used as biological evidence.

The `expected` directory contains seven scientific CSV tables generated
with the reference R backend after independent count assertions in the
[fixture generator](https://github.com/castrattonDSU/FibrTraQ/blob/main/data-raw/generate-known-answer-example.R).
`ftq_run_example()` compares new outputs against those tables and resolves
all 18 declarations in `known-answer-claims.csv` from both the stored and
newly generated tables. Backend identity and installation-path provenance
are normalized before comparison; all scientific values, labels,
denominators, and row ordering must match exactly.

The fixture was authored for FibrTraQ and is distributed under the same
GNU General Public License, version 2 or later, as the package.
