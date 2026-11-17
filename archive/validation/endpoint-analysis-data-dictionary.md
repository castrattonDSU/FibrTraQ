# Endpoint analysis data dictionary

FibrTraQ's primary per-voxel and ROI target counts use the selected
target-assignment mode: full streamline, endpoint voxels, or local search.
The additive endpoint-audit tables deliberately retain full-path target
intersection as their denominator so endpoint membership can be interpreted
against the same whole-path estimand in every mode. These quantities must not
be substituted for one another in tables, figures, or manuscript text.

For endpoint-voxel and local-search runs, the ROI-summary CSV persists both:

- `unique_target_streamlines:<target>`: the primary count under the selected
  assignment mode.
- `unique_full_path_target_streamlines:<target>`: the independently retained
  whole-path count used to validate the endpoint-audit denominator.

For full-streamline runs, `unique_target_streamlines:<target>` already is the
whole-path count, so the additional `unique_full_path_*` field is omitted to
preserve the established output schema.

## Endpoint summary

`*_fibrtraq_endpoint-summary.csv` has one row per target and endpoint tolerance.

- `tolerance`: `exact` is primary ROI membership. `dilated1` is a Chebyshev
  radius-one, 26-neighbor voxel-grid sensitivity analysis and is not a uniform
  millimeter tolerance.
- `unique_focal_intersecting_streamlines`: unique streamlines intersecting the
  focal ROI under the selected FibrTraQ voxelization.
- `target_intersection_streamlines`: the established full-path target-positive
  count among focal-intersecting streamlines.
- `focal_endpoint_streamlines`: focal-intersecting streamlines with at least one
  endpoint in the focal ROI at the stated tolerance.
- `paired_endpoint_streamlines`: all focal-positive streamlines with opposite
  focal and target endpoint memberships at the stated tolerance. This raw pair
  quantity is monotonic from exact to one-voxel membership and can include a
  streamline whose endpoint also has shared focal-target membership.
- `target_intersection_with_paired_endpoints`: full-path target intersections
  that also satisfy the raw paired-endpoint definition.
- `shared_mask_ambiguous`, `paired_focal_target_endpoints`,
  `target_endpoint_only`, `focal_endpoint_only`, and
  `neither_relevant_endpoint`: mutually exclusive classes that partition
  `target_intersection_streamlines` in that priority order.
- `paired_focal_target_endpoints` is therefore the unambiguous exclusive class.
  It can be smaller than `target_intersection_with_paired_endpoints`; any
  difference is contained in `shared_mask_ambiguous`.
- Percentage column names state their denominator. A missing percentage means
  its denominator was zero; it must not be converted to zero.
- `focal_target_overlap_voxels` quantifies ROI overlap at each tolerance. Shared
  endpoint membership is reported as ambiguity rather than forced into a pair.

## Penetration

`*_fibrtraq_endpoint-penetration.csv` has one row per target and observed number
of unique exact target voxels visited. `streamline_count` values sum to the
target's full-path intersection count. Under `sample_points`, penetration counts
only voxels represented by stored tractography points and can undercount a
continuous path. Under `segment_interpolation`, the selected interpolated voxel
path is used.

## Endpoint context

`*_fibrtraq_endpoint-context.csv` always contains `reference_grid` and contains
one additional row for every explicitly supplied binary anatomical context.

- `endpoint_count` counts one endpoint for one-point streamlines and two for all
  other valid streamlines.
- `endpoints_inside` and `endpoints_outside` sum to `endpoint_count`.
- Reference-grid rows report axis-specific boundary crossings and distinguish
  exactly one voxel from more than one voxel outside the grid.
- Binary-mask rows report endpoint and streamline containment but do not apply
  an automatic biological pass threshold.
- Intensity images are never thresholded to create endpoint contexts.

## Validation status

`validation$passed` covers computational invariants and reconciliation of
endpoint-audit denominators with the applicable persisted whole-path target
counts. `biological_qc_status` is separate:

- `pass`: an anatomical context was supplied and no severe grid/non-finite or
  shared-mask ambiguity signal was detected.
- `review`: a severe grid/non-finite or shared-mask ambiguity signal requires
  interpretation.
- `unavailable`: no explicit anatomical context was supplied or endpoint tables
  were unavailable.

Endpoint analysis strengthens biological interpretation but does not establish
tractogram/reference provenance. The reconstruction reference geometry and the
provenance of any brain/cord mask remain external evidence requirements.
