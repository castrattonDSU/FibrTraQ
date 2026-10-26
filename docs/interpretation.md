# Interpreting FibrTraQ results

[Back to the quick start](../README.md) · [Advanced workflows](advanced-workflows.md) · [Complete output dictionary](../inst/output-guide/OUTPUTS_README.md)

## RSC measures focal-voxel coverage

For one target:

```text
RSC (%) = 100 × focal voxels with at least one association to this target
               / all focal voxels
```

The denominator includes zero-assignment focal voxels. A focal voxel can have
associations with multiple targets, so target RSC values can overlap and do not
have to sum to 100%. The `local` connection row can overlap target-positive
rows as well. An RSC of 100% means every focal voxel has an association with
that target under the selected method; it does not mean every streamline is
assigned to that target or that the association is biologically valid.

Under `first_intersection`, a zero-assignment voxel can still be traversed by a
streamline allocated to another focal voxel. Read RSC and assignment coverage
with the recorded focal-allocation method.

## Per-voxel counts are not ROI-level unique streamlines

Under the recommended default,
`focal_assignment_mode = "all_focal_voxels"`, one streamline can occupy
multiple focal voxels. Therefore:

- a per-voxel `tracks` value is a focal-voxel/streamline association count;
- summing per-voxel values does not produce a unique-streamline total; and
- unique focal and target-positive streamline totals should be read from the
  ROI-summary and target-summary outputs.

`ftq_target_summary_tables()` keeps unique-streamline and voxel-association
estimands in separate, explicitly labeled tables.

## Target-assignment modes

The direct-call `assignment_mode` argument (stored in manifests and provenance
as `target_intersection_mode`) controls the primary target columns,
target-positive unique-streamline counts, and target-dependent fate
classification:

| Canonical value | Primary target-assignment rule |
|---|---|
| `"full_streamline"` | Backward-compatible all-voxel mode: any selected path voxel in a target records that target connection |
| `"endpoint_voxels"` | The two endpoint voxels, under the first focal-positive method in the ordered `discretization` candidates, must resolve unambiguously to the focal ROI and one target ROI |
| `"local_search"` | Each precise endpoint is assigned independently to the nearest labelled focal/target voxel centre within `search_radius_mm`, after which an unambiguous focal-target pair is required |

For `"local_search"`, distances are scanner-space Euclidean distances in
millimetres from the continuous endpoint coordinate to voxel centres, and the
radius is inclusive. The default is 2 mm. This mode therefore requires
trustworthy voxel-to-scanner geometry; direct `FibrTraQ()` callers using plain
arrays must supply `ref_image`. Dimensions are checked, but the caller must
establish that it is the exact `fNode` grid because the focal array has no
affine metadata.
NIfTI metre and micrometre spatial units are converted to millimetres. When
the header reports unknown spatial units, FibrTraQ assumes millimetres, emits
a warning, and records that assumption in preflight and result provenance. If
an endpoint's nearest distances to distinct labels are equal within the
implementation's relative squared-distance tolerance
(`1e-12 * max(1, a, b)`), or a candidate voxel belongs to overlapping
focal/target labels, the endpoint is ambiguous and no primary connection is
recorded. Multiple equally near voxels carrying the same label remain one
unambiguous label assignment.

For a local-search run, FibrTraQ also writes a dedicated
`*_fibrtraq_local-search-assignments.csv` audit. Use
`ftq_local_search_tables(result)` to retrieve its streamline-level rows and
summaries in memory. The audit partitions each assigned cNode connection into:

- **also intersects cNode**, where the selected voxelized path enters the
  assigned target mask; and
- **proximity only**, where both ends satisfy the local-search rule but the
  selected path does not enter the assigned target mask.

These are different geometric statements. A proximity-only assignment is a
local-search connection, but it must not be described as a direct target-mask
intersection. The audit also reports end-voxel distance quantiles so the
effective tolerance can be evaluated instead of inferred from the configured
maximum radius alone.

Endpoint-voxel and local-search assignments are joint over the focal mask and
the complete supplied target set. Adding, removing, or changing a target can
therefore change another target's assignment by introducing a nearer label, an
equal-distance tie, or an overlapping label. Preserve the complete target list
in provenance and method comparisons.

Paper-aligned aliases are normalized in saved provenance:
`"all_voxel"` and `"all_voxels"` map to `"full_streamline"`; `"end_voxel"`,
`"end_voxels"`, and `"endpoint_voxel"` map to `"endpoint_voxels"`. Prefer the
canonical values in new manifests.

These choices correspond to the broad assignment strategies compared by
[Yeh et al. (2019)](https://doi.org/10.1016/j.neuroimage.2019.05.005).
The 2 mm default follows the local-search radius reported in that study; it is
not assumed to be universally optimal.
FibrTraQ's local-search rule above is an explicit, reproducible operational
contract for details not fully specified in that paper; it should not be
described as an exact reproduction of the paper's RC3 implementation.

Changing target assignment does not change which streamlines qualify through
the focal path or how qualifying streamlines are allocated to focal voxels.
Those established semantics remain centralized in `FibrTraQ()`, including
`focal_assignment_mode` and `streamline_voxelization_mode`. The default
`"full_streamline"` preserves established counting behavior and the primary
result-table schema.

Primary-negative fate counters currently combine several reasons for no
assigned target, including no label within range, a distinct-label tie,
overlapping labels, or endpoints that do not form a focal-target pair. Preserve
the complete target set and inspect spatial QC when those reasons matter.

The output field and fate named `local` are unrelated to `"local_search"`.
`local` is a non-exclusive classification based on whether a focal-positive
streamline stays inside `localROI`; `"local_search"` is a target-assignment
method based only on the two endpoints.

## Publication report summaries use unique streamlines

The standard HTML report adds two compact, publication-facing tables without
changing any established result schema:

- The connection summary contains one `local` row and one row per cNode. It
  reports unique streamline connections, connected fNode voxels, RSC, and the
  total number of fNode voxels.
- The report fate summary partitions each unique focal-intersecting streamline
  exactly once into `local`, `leave_with_target`, or `leave_no_target`.

The unique-streamline counts in these report tables come from ROI-summary
counters and must not be reconstructed by summing per-voxel values. The
established RSC CSV retains association totals and the same focal-voxel RSC
measure; its association columns are not unique counts. The established fate
CSV uses voxel-weighted associations, while the report fate summary uses
unique streamlines when their ROI-level provenance is available. The output
dictionary identifies the explicit association-based fallback for standalone
reports that lack that provenance.

Density heatmaps in the HTML report use a display-only 99th-percentile cap so
isolated high-density voxels do not hide the remaining spatial structure. Each
connection-density row has its own labeled color scale. Consequently, color
intensity is comparable within a row but not between differently colored rows.
The underlying CSV and NIfTI values are never clipped or transformed by this
display setting.

## The additive endpoint audit answers a separate question

- **Full-path intersection audit:** Did any selected path voxel of a
  focal-positive streamline occupy the target mask?
- **Exact paired endpoints:** Did opposite streamline endpoints lie
  unambiguously in the focal and target masks?
- **One-voxel paired endpoints:** Does that result change after a fixed
  26-neighbor, one-voxel dilation?
- **Penetration:** How many unique exact target voxels did a full-path
  intersection-positive streamline visit?

The existing `endpoint_analysis = "exact_and_dilated1"` tables remain additive
QC outputs. They retain the full-path intersection denominator even when the
primary saved `target_intersection_mode` is `"endpoint_voxels"` or
`"local_search"`;
they do not redefine the selected primary assignment. The one-voxel result is
a grid sensitivity analysis, not a uniform physical-distance tolerance and not
the local search.

## Computational and biological gates are separate

`validation_passed` means the package's computational reconciliation checks
passed. `biological_qc_status` may still be `review` or `unavailable` because
tractogram provenance, anatomical containment, or spatial alignment requires
human or external evidence. Do not convert computational success into a
biological claim without that evidence.

## The 3-D viewer is a QC sample

The interactive viewer preserves a target-balanced quota subset for visual
review. For each cNode it accepts positive streamlines in tractogram scan order
until that cNode's quota or the configured scan limit is reached. It is not a
random sample, is not prevalence weighted, is not balanced by end-voxel class,
and does not display every streamline. Quantitative claims must come from the
full-run tables and maps.

The initial camera uses an anatomical RAS display basis with superior up and a
rotating R/L, A/P, S/I compass. End voxels are shown by default. Dense fNode and
cNode hit markers, path-segment modes, and end-voxel class filters remain
available under advanced controls so they do not obscure the initial geometry
review.
