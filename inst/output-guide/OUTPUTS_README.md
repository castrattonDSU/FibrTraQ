# FibrTraQ output guide

Definitions for the installed package's public file workflows. Keep this guide with the complete results folder.

## Start with Results

Open the study or single-run report index and choose Results. Read the subject, focal region, assignment method, computational status and biological QC separately. The RSC summary, core maps and complete voxel table are the main results; advanced tables remain available for audit.

The complete original *_fibrtraq_per-voxel-hits.csv contains every focal voxel. The readable *_all-voxels.csv adds explicit array coordinates; use *_all-voxels-columns.csv to identify coordinate and target roles safely. Ranked CSVs contain every focal voxel, including zero-count voxels: descending selected metric, then increasing x, y, z for ties. Historical top-50 files contain only a preview.

Use the original per-voxel-hits CSV as the general analysis input. The all-voxel file is a readable export: a plain read.csv() retains coordinate columns, which must not be passed to target-inference functions as if they were target counts. FibrTraQ file readers recognize the named all-voxel CSV together with its companion column-role record and restore the canonical fields.

Copy the complete result folder and all subfolders, then extract it before opening HTML. HTML uses portable relative links; recorded source paths in CSV/JSON preserve provenance from the generating workstation. JSON is machine-readable text, not PDF; CSV is a table, not an Excel workbook.

## Report locations

The main Results report is under reports/<run>/; its images and embedded viewer are in report-figures/. Complete ranked CSVs are in tables/ beside the report. The supporting Target report and its five PNG figures are in the focal-region figures/ folder. The primary-target lookup CSV is in the focal-region folder. Standard runs do not create exploratory-report/ or context/ folders.

The Target report retains its calculated summaries in memory and displays the method and denominators. For further research, ftq_exploratory_summary() can explicitly export the detailed derived tables and maps to a separate directory. Existing research exports and historical result folders remain readable.

## Counting populations and missing values

A unique streamline count counts a tractogram record once within its stated population. A focal-voxel/streamline association counts its allocation to a particular focal voxel; summing across voxels can count a streamline repeatedly. Target-association totals can also repeat across targets. Never infer biological sample size from voxels, streamlines or workflow runs.

RSC (%) = 100 * focal voxels with a positive count for this target / ALL focal voxels, including zero-assignment voxels. It is not a percentage of streamlines. A focal voxel can contribute to several targets, even when each streamline is assigned to only one target, so RSC percentages do not need to sum to 100%. Local means every selected path voxel lies in the configured local ROI, which may differ from the focal ROI. Local and target-positive connection rows can overlap. The three fate rows are exclusive because local takes precedence, and the two nonlocal rows split by target assignment.

All focal voxels are retained. Under all_focal_voxels, a streamline is allocated to each distinct occupied focal voxel. Under first_intersection, only its first selected focal voxel receives the count: a zero-assignment voxel can still be traversed. Do not equate zero-assignment coverage with biological absence or endpoint coverage.

NA means unavailable, inapplicable or unevaluable as specified for each field, not zero. Several legacy fraction fields store 0 for zero denominators: the dictionary identifies them. Exclude unevaluable voxels when interpreting distributions. The additive local-fraction CSV uses NA for zero tracks; its NIfTI stores 0 with a separate validity mask.

## Method and endpoint audit

full_streamline assigns targets from selected path voxels; endpoint_voxels requires an unambiguous focal/target endpoint pair; local_search requires an unambiguous nearest-label pair within the inclusive physical radius. The complete target set affects endpoint/local-search assignment through overlap, ties and nearer labels. Preserve that set and the recorded method when rebuilding reports.

The additive endpoint-summary and penetration tables always audit exact full-path target intersection. Exact endpoint membership and dilated1 (Chebyshev radius one voxel, including 26 neighbours) are separate from physical-radius local search. Broad paired_endpoint counts may include shared-mask ambiguity; the unambiguous paired_focal_target_endpoints class excludes it. Use the numerator named in each percentage definition.

Endpoint coverage counts all valid mapped streamlines before focal filtering. 100% inside the image grid means endpoints lie inside image bounds; it does not establish anatomy, continuous path coverage, local classification or target assignment. With trusted masks, anatomical endpoint containment remains QC evidence rather than proof of biological validity.

Penetration is the number of unique exact target voxels occupied by the selected full path. It is not physical depth. Histogram bins are sparse; their counts sum to full-path target-intersection streamlines. Under sample_points, sparse stored points can miss intervening voxels. The five exclusive endpoint classes likewise partition that same target-intersection population.

## Maps, coordinates and image lookup

fNode_Coords and readable x/y/z values are one-based image-array indices, not scanner millimetres. Image affines and declared spatial units supply the physical transformation. Existing target names, including x, y or z, are preserved; the companion column-role table identifies renamed coordinate fields.

Core maps in maps/ retain the focal reference grid and affine: tracks (assigned streamline count), tracks-binary (1 if tracks > 0), local, leave-no-target, leave-with-target, and target-{label}. Outside the focal mask these maps are zero. The original focal mask is necessary to distinguish outside-mask zeros from zero-assignment focal voxels. Display scaling affects previews only.

The additive local-fraction map is local / tracks on positive-assignment focal voxels, with stored zero elsewhere. The accompanying local-fraction-valid map is 1 only where tracks > 0, including zero-fraction but evaluable voxels. Read them together. These maps use the configured local ROI, not an assumed focal-only region.

Exploratory maps are available only through the explicit research export ftq_exploratory_summary(..., output_dir = ..., template_image = ..., write_maps = TRUE). Standard reports do not write these maps or their supporting CSVs. Research maps retain the corresponding exploratory-per-voxel quantities: primary_target_index, primary_target_hits, primary_target_fraction, target_diversity, target_entropy_bits, target_entropy_fraction, total_target_hits, target_hits_per_assigned_track, leave_no_target_fraction and leave_with_target_fraction. They are supporting views. Category index 0 means no target; positive indices use the primary-target LUT. Ties choose the first target column. A stored zero entropy can mean no data or a concentrated distribution, so inspect target/assignment coverage.

Tractogram density context is whole-tractogram streamline occupancy per reference voxel (at most once per streamline per voxel). It is a context image, not MRI, a brain mask or anatomical tissue volume. Core maps and context maps have different counting populations. The density image is saved with the core maps; its summary and settings are retained in the run metadata. PNG figures are previews; their portable download links provide full-size images, while NIfTI files preserve the full voxel values.

The 3-D viewer is a quota-limited, scan-order subset selected by target assignment. It is neither random nor prevalence weighted. Its paths and endpoint markers support spatial inspection; they must not be used to estimate target prevalence or whole-tractogram completeness.

## Methods, validation and supporting records

*_metadata.json preserves generation time, loaded package identity, normalized manifest, input paths/checksums, compact preflight, counting summaries, endpoint method, local-search geometry and output paths. *_validation.json preserves computational check results, reconciliation and expected-baseline checks where requested, plus separate biological QC warnings/status. Readable Methods and Validation sections expose these records without replacing their JSON source.

Prepared input, union and selected-volume CSVs document what was used or derived. The prepared/ folder is created only when an image or mask derivative is needed. Optional unavailable inputs are not analysis failures. Run indexes identify every subject/focal analysis and its outputs; record and inspect analysis versus presentation status separately. Local output paths and subject labels may contain information that needs review before sharing publicly.

The bundled synthetic example compares seven scientific tables and 18 independently declared claims exactly. Four HTML artifacts are checked for existence and nonzero size. These checks do not compare rendered HTML exactly, exercise browser interaction or establish biological validity.

The families below describe standard file-workflow outputs, explicitly requested research exports, historical fields and the public grid/validation diagnostic CSVs. Target labels and added coordinates are explicit parameterized roles, not inferred numeric columns. Optional columns/whole tables appear only when their source feature is enabled. All legacy aliases remain documented so saved analyses remain readable.

## per_voxel

Files or returned table: `*_fibrtraq_per-voxel-hits.csv`.

### `fNode_Coords`

Coordinate key x_y_z, using one-based focal-image array indices.

- Units: voxel index key.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Every focal voxel has a key, including zero-assignment voxels.
- Method: selected assignment mode.
- Alias of: none.

### `tracks`

Streamlines allocated to this focal voxel; each streamline contributes at most once to a voxel.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `local`

Allocated streamlines whose entire selected voxelized path is inside the configured local ROI. Local status can overlap target-positive status.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target`

Allocated streamlines that are nonlocal and receive no target under the selected assignment rule.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `{target}`

Allocated streamlines assigned to the named target. The target name is supplied by the user; counts across target columns may overlap.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `none`

1 when tracks is zero, otherwise 0. This is zero assignment, not necessarily absence of a traversing streamline under first_intersection.

- Units: binary indicator.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means at least one assigned streamline; 1 means none assigned.
- Method: selected assignment mode.
- Alias of: none.

### `cNode_Coords`

Optional legacy coordinate annotation; not a target-count column.

- Units: coordinate annotation.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: legacy input only.
- Alias of: none.

## all_voxels

Files or returned table: `tables/*_all-voxels.csv`.

### `fNode_Coords`

Coordinate key x_y_z, using one-based focal-image array indices.

- Units: voxel index key.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Every focal voxel has a key, including zero-assignment voxels.
- Method: selected assignment mode.
- Alias of: none.

### `tracks`

Streamlines allocated to this focal voxel; each streamline contributes at most once to a voxel.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `local`

Allocated streamlines whose entire selected voxelized path is inside the configured local ROI. Local status can overlap target-positive status.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target`

Allocated streamlines that are nonlocal and receive no target under the selected assignment rule.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `{target}`

Allocated streamlines assigned to the named target. The target name is supplied by the user; counts across target columns may overlap.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `none`

1 when tracks is zero, otherwise 0. This is zero assignment, not necessarily absence of a traversing streamline under first_intersection.

- Units: binary indicator.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means at least one assigned streamline; 1 means none assigned.
- Method: selected assignment mode.
- Alias of: none.

### `cNode_Coords`

Optional legacy coordinate annotation; not a target-count column.

- Units: coordinate annotation.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: legacy input only.
- Alias of: none.

### `{coordinate}`

Added x/y/z array index. Consult the companion all-voxels-columns.csv: existing target names are retained; conflicting coordinates use voxel_x, voxel_y, voxel_z and a unique suffix.

- Units: one-based voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; never infer the role from a column name.
- Method: focal image grid.
- Alias of: none.

## all_voxels_columns

Files or returned table: `tables/*_all-voxels-columns.csv`.

### `column`

Exact column name in the companion all-voxels.csv.

- Units: text.
- Population: one exported column.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive readable export.
- Alias of: none.

### `role`

Explicit role: coordinate, identifier, count or target. Use this to distinguish coordinates from targets named x, y or z.

- Units: text.
- Population: one exported column.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive readable export.
- Alias of: none.

### `source_column`

Source field in the original per-voxel table; coordinates are derived from fNode_Coords.

- Units: text.
- Population: one exported column.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive readable export.
- Alias of: none.

### `units`

Measurement or indexing unit for this column.

- Units: text.
- Population: one exported column.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive readable export.
- Alias of: none.

### `description`

Human-readable definition of this column.

- Units: text.
- Population: one exported column.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive readable export.
- Alias of: none.

## local_fraction

Files or returned table: `tables/*_local-fraction.csv`.

### `fNode_Coords`

Coordinate key x_y_z, using one-based focal-image array indices.

- Units: voxel index key.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Every focal voxel has a key, including zero-assignment voxels.
- Method: selected assignment mode.
- Alias of: none.

### `local_fraction`

local / tracks.

- Units: fraction.
- Population: one focal voxel.
- Denominator: tracks.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: selected assignment mode.
- Alias of: none.

### `has_assigned_streamlines`

TRUE when tracks > 0; identifies where local_fraction is evaluable.

- Units: logical.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: FALSE means zero assignment and unevaluable fraction.
- Method: selected assignment mode.
- Alias of: none.

## rsc

Files or returned table: `*_fibrtraq_rsc-summary.csv`.

### `node`

Target label; the RSC source CSV contains target rows.

- Units: text.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `target_associations`

Sum of this target's per-voxel counts; one streamline may contribute at several focal voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxel_streamline_target_associations`

Same count as target_associations.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: target_associations.

### `tracks`

Legacy alias of target_associations; not a unique streamline total.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: target_associations.

### `voxels_connected`

Focal voxels with a positive count for this target.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `total_focal_voxels`

All focal voxels, including those with zero assigned streamlines.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `total_origin_voxels`

Legacy alias of total_focal_voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: total_focal_voxels.

### `rsc`

voxels_connected / total_focal_voxels.

- Units: fraction.
- Population: one target.
- Denominator: total_focal_voxels.
- Missing/zero: 0 means no numerator events when evaluable; The legacy empty-table case stores 0; it is unevaluable.
- Method: selected assignment mode.
- Alias of: none.

### `rsc_percent`

100 * voxels_connected / total_focal_voxels.

- Units: percent.
- Population: one target.
- Denominator: total_focal_voxels.
- Missing/zero: 0 means no numerator events when evaluable; The legacy empty-table case stores 0; it is unevaluable.
- Method: selected assignment mode.
- Alias of: none.

## fate

Files or returned table: `*_fibrtraq_fate-summary.csv (voxel-weighted)`.

### `fate`

Exclusive voxel-weighted class: local, leave_no_target or leave_with_target.

- Units: text.
- Population: one fate class.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `tracks`

Sum of per-voxel counts for this fate; leave_with_target = tracks - local - leave_no_target.

- Units: count.
- Population: focal-voxel/streamline associations.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `percent`

100 * tracks in this row / sum of tracks across the three fate rows.

- Units: percent.
- Population: focal-voxel/streamline associations.
- Denominator: sum of tracks across the three fate rows.
- Missing/zero: 0 means no numerator events when evaluable; Legacy zero-denominator storage is 0, an unevaluable value.
- Method: selected assignment mode.
- Alias of: none.

## roi_summary

Files or returned table: `*_fibrtraq_roi-summary.csv`.

### `total_input_streamlines`

Tractogram's declared total, if known; it can differ from the number processed when a limit is used.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `processed_streamlines`

Streamlines processed during the tractogram pass, before focal-positive filtering.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_focal_intersecting_streamlines`

Unique streamlines with at least one selected focal voxel.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_local_streamlines`

Focal-positive streamlines entirely local, irrespective of target assignment.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_leave_no_target_streamlines`

Focal-positive nonlocal streamlines without a target assignment.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_leave_with_target_streamlines`

Focal-positive nonlocal streamlines with one or more assigned targets.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_focal_intersecting_with_target_streamlines`

Focal-positive streamlines with any target assignment, including local ones.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_local_target_positive_streamlines`

Local and target-positive focal streamlines: one of four exclusive local-by-target groups.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_local_target_negative_streamlines`

Local and target-negative focal streamlines: one of four exclusive local-by-target groups.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_nonlocal_target_positive_streamlines`

Nonlocal and target-positive focal streamlines; same as unique_leave_with_target_streamlines.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: unique_leave_with_target_streamlines.

### `unique_nonlocal_target_negative_streamlines`

Nonlocal and target-negative focal streamlines; same as unique_leave_no_target_streamlines.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: unique_leave_no_target_streamlines.

### `unique_any_target_positive_streamlines`

Focal-positive streamlines with any assigned target; same as unique_focal_intersecting_with_target_streamlines.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: unique_focal_intersecting_with_target_streamlines.

### `unique_any_target_negative_streamlines`

Focal-positive streamlines without any assigned target, combining local and nonlocal groups.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_exclusive_fate_sum`

Sum of the four exclusive local-by-target groups; should equal unique_focal_intersecting_streamlines.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `total_focal_voxel_streamline_associations`

Sum of per-voxel tracks, with possible repeat contributions by one streamline across focal voxels.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `total_focal_voxel_streamline_target_associations`

Sum over all focal voxels and targets; a streamline can contribute across both voxels and targets.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_target_streamlines:{target}`

Unique focal-positive streamlines assigned to this target under the primary method.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_full_path_target_streamlines:{target}`

Unique focal-positive streamlines whose selected full path enters this exact target mask; additive comparator in endpoint/local-search runs.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_exclusive_fate_sum_matches_focal`

Arithmetic equality of unique_exclusive_fate_sum and unique_focal_intersecting_streamlines.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_assignment_mode`

all_focal_voxels counts each occupied focal voxel; first_intersection allocates a streamline only to its first selected focal voxel.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `streamline_voxelization_mode`

sample_points uses stored tractogram points; segment_interpolation densifies consecutive segments before voxelization.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `target_intersection_mode`

Primary assignment: full_streamline, endpoint_voxels or local_search. These change target assignment, not the separate full-path endpoint audit.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `local_roi_defaulted_to_focal`

Whether the focal mask itself defined the local ROI; FALSE means a separately supplied or prepared local region was used.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `selected_track_space`

Coordinate mapping actually used for tractogram points.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `track_voxel_size_mm`

Comma-separated TrackVis voxel sizes used for coordinate conversion, in millimetres.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `selected_axis_flip`

Comma-separated logical axis flips in x,y,z order.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `selected_discretization`

Ordered voxel discretization candidates. For any, a focal-positive path uses the first candidate that qualifies it.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `reverse_streamlines`

Whether stored point order was reversed before counting.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `backend`

Counting implementation used: R or compiled backend; not a biological method.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `mapping_selection_rule`

Rule selecting automatic coordinate mapping; target/local hits do not select or break ties.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `local_search_spatial_units`

Spatial-unit declaration in the reference image header.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `local_search_spatial_units_interpretation`

Recorded interpretation of those units, including an explicit millimetre assumption for unknown units.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `summary_available`

FALSE in a fallback record when a reloaded per-voxel CSV lacks unique-streamline provenance.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `note`

Explanation of unavailable summary provenance.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `local_search_radius_mm`

Inclusive maximum endpoint-to-labelled-voxel-centre distance.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search only.
- Alias of: none.

### `local_search_spatial_units_scale_to_mm`

Multiplicative factor converting declared image-space distance to millimetres.

- Units: conversion factor.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search only.
- Alias of: none.

## connection_summary

Files or returned table: `tables/*_connection-summary.csv`.

### `cNode`

Target label, or local for the configured local region. Local can overlap target-positive rows.

- Units: text.
- Population: one target/local row.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_streamline_connections`

ROI-level unique streamline count for this target or local class; NA if unavailable from saved provenance.

- Units: count.
- Population: one target/local row.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `connected_fNode_voxels`

Focal voxels with a positive target/local count.

- Units: count.
- Population: one target/local row.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `total_fNode_voxels`

All focal voxels, including zero-assignment voxels.

- Units: count.
- Population: one target/local row.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `rsc_percent`

100 * connected_fNode_voxels / total_fNode_voxels.

- Units: percent.
- Population: one target/local row.
- Denominator: total_fNode_voxels.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: selected assignment mode.
- Alias of: none.

## report_fate

Files or returned table: `tables/*_report-fate-summary.csv`.

### `fate`

Exclusive class: local, leave_no_target or leave_with_target. Read summary_basis before interpreting counts.

- Units: text.
- Population: one fate class.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `count`

Count for this class in the stated summary_basis.

- Units: count.
- Population: summary_basis.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `denominator_count`

All focal-positive unique streamlines, or all focal-voxel/streamline associations in the explicit fallback.

- Units: count.
- Population: summary_basis.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `summary_basis`

unique_focal_intersecting_streamlines when available; focal_voxel_streamline_associations is a clearly declared fallback.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `percent_of_denominator`

100 * count / denominator_count.

- Units: percent.
- Population: summary_basis.
- Denominator: denominator_count.
- Missing/zero: 0 means no numerator events when evaluable; Legacy zero-denominator storage is 0; the fraction is unevaluable.
- Method: selected assignment mode.
- Alias of: none.

## top_voxels

Files or returned table: `tables/*_ranked_*_voxels.csv (all focal voxels; historical top-N previews remain readable)`.

### `tracks`

Streamlines allocated to this focal voxel; each streamline contributes at most once to a voxel.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `local`

Allocated streamlines whose entire selected voxelized path is inside the configured local ROI. Local status can overlap target-positive status.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target`

Allocated streamlines that are nonlocal and receive no target under the selected assignment rule.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `{target}`

Allocated streamlines assigned to the named target. The target name is supplied by the user; counts across target columns may overlap.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `x`

One-based focal-image array index along x; not a scanner-space distance.

- Units: voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; zero and NA are invalid.
- Method: focal image grid.
- Alias of: none.

### `y`

One-based focal-image array index along y; not a scanner-space distance.

- Units: voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; zero and NA are invalid.
- Method: focal image grid.
- Alias of: none.

### `z`

One-based focal-image array index along z; not a scanner-space distance.

- Units: voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; zero and NA are invalid.
- Method: focal image grid.
- Alias of: none.

### `leave_with_target`

tracks - local - leave_no_target; nonlocal allocated streamlines with any assigned target.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `{derived_leave_with_target}`

Same derived nonlocal target-positive count as tracks - local - leave_no_target. If a genuine target is named leave_with_target, the derived field uses derived_leave_with_target, repeatedly adding derived_ until unique; preserve the target's original label.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target.1`

Duplicate leave_no_target column in older top-50 previews; retained in historical files, omitted from newly generated ranked exports.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: historical preview only.
- Alias of: leave_no_target.

## endpoint_summary

Files or returned table: `*_fibrtraq_endpoint-summary.csv`.

### `target`

Target whose exact full-path intersection population is audited.

- Units: text.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `tolerance`

exact or dilated1; the latter includes the 26 neighbouring voxels and is not a uniform millimetre radius.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `tolerance_definition`

Human-readable definition of the endpoint-membership tolerance.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `unique_focal_intersecting_streamlines`

All unique focal-positive streamlines in the run.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_streamlines`

Unique focal-positive streamlines whose selected full path occupies at least one exact target voxel.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_streamlines`

Focal-positive streamlines with at least one focal endpoint at this tolerance, irrespective of target intersection.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_endpoint_streamlines`

Focal-positive streamlines with opposite focal and target endpoints, including ambiguous shared-mask membership and pairs without an exact target intersection.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_with_paired_endpoints`

Exact target-intersection streamlines with opposite focal and target endpoints; this broad count can include shared-mask ambiguity.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_with_any_target_endpoint`

Exact target-intersection streamlines with at least one target endpoint, including ambiguous cases.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_with_any_focal_endpoint`

Exact target-intersection streamlines with any focal endpoint: shared-mask-ambiguous, paired and focal-only classes combined.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_with_target_intersection_streamlines`

Focal-positive streamlines with a focal endpoint and an exact target intersection, including ambiguous membership.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `shared_mask_ambiguous`

Target-intersection streamlines with an endpoint in both focal and target masks at this tolerance; precedence class in the five-class partition.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_focal_target_endpoints`

Unambiguous opposite focal/target endpoints among exact target-intersection streamlines; excludes shared-mask ambiguity.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_endpoint_only`

Unambiguous target endpoint present, no focal endpoint, among target-intersection streamlines.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_only`

Unambiguous focal endpoint present, no target endpoint, among target-intersection streamlines.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `neither_relevant_endpoint`

Neither relevant endpoint in the five-class target-intersection partition.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_mask_voxels`

Focal-mask voxel count at the displayed exact/dilated tolerance.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_mask_voxels`

Target-mask voxel count at the displayed exact/dilated tolerance.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_target_overlap_voxels`

Voxel overlap of focal and target masks at this tolerance.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_percent_of_unique_focal_intersecting_streamlines`

100 * target_intersection_streamlines / unique_focal_intersecting_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: unique_focal_intersecting_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_percent_of_target_intersections`

100 * target_intersection_with_paired_endpoints / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_percent_of_focal_endpoint_streamlines`

100 * paired_endpoint_streamlines / focal_endpoint_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: focal_endpoint_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_endpoint_percent_of_target_intersections`

100 * target_intersection_with_any_target_endpoint / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_percent_of_target_intersections`

100 * target_intersection_with_any_focal_endpoint / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_percent_of_focal_endpoint_streamlines`

100 * focal_endpoint_with_target_intersection_streamlines / focal_endpoint_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: focal_endpoint_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_target_overlap_percent_of_smaller_mask`

100 * focal_target_overlap_voxels / min(focal_mask_voxels, target_mask_voxels).

- Units: percent.
- Population: one target/tolerance.
- Denominator: min(focal_mask_voxels, target_mask_voxels).
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `shared_mask_ambiguous_percent_of_target_intersections`

100 * shared_mask_ambiguous / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_endpoint_only_percent_of_target_intersections`

100 * target_endpoint_only / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_only_percent_of_target_intersections`

100 * focal_endpoint_only / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `neither_relevant_endpoint_percent_of_target_intersections`

100 * neither_relevant_endpoint / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `voxel_size_x_mm`

Reference voxel spacing along x after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `voxel_size_y_mm`

Reference voxel spacing along y after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `voxel_size_z_mm`

Reference voxel spacing along z after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `streamline_voxelization_mode`

sample_points uses stored tractogram points; segment_interpolation densifies consecutive segments before voxelization.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_mode`

Always full_streamline in this additive audit, even when the run's primary assignment uses endpoints or local search.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

## endpoint_penetration

Files or returned table: `*_fibrtraq_endpoint-penetration.csv`.

### `target`

Target whose exact full-path intersection population is audited.

- Units: text.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_streamlines`

Unique focal-positive streamlines whose selected full path occupies at least one exact target voxel.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `streamline_voxelization_mode`

sample_points uses stored tractogram points; segment_interpolation densifies consecutive segments before voxelization.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `unique_target_voxels`

Exact distinct target voxels visited by a selected full path; histogram bin value, not depth in mm.

- Units: count.
- Population: one target/penetration bin.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `streamline_count`

Unique target-intersection streamlines in this penetration bin.

- Units: count.
- Population: one target/penetration bin.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `percent_of_target_intersections`

100 * streamline_count / target_intersection_streamlines.

- Units: percent.
- Population: one target/penetration bin.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `penetration_interpretation`

Stored-point or segment-interpolated penetration definition; sparse stored points may miss intervening target voxels.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

## penetration_summary

Files or returned table: `tables/*_endpoint-penetration-summary.csv`.

### `target`

Target whose exact full-path intersection population is audited.

- Units: text.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_streamlines`

Unique focal-positive streamlines whose selected full path occupies at least one exact target voxel.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `one_voxel_touch_streamlines`

Target-intersection streamlines visiting exactly 1 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `one_voxel_touch_percent_of_target_intersections`

100 * one_voxel_touch_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `two_voxel_only_streamlines`

Target-intersection streamlines visiting exactly 2 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `two_voxel_only_percent_of_target_intersections`

100 * two_voxel_only_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `three_to_four_voxels_streamlines`

Target-intersection streamlines visiting 3 or 4 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `three_to_four_voxels_percent_of_target_intersections`

100 * three_to_four_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `five_to_nine_voxels_streamlines`

Target-intersection streamlines visiting 5 through 9 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `five_to_nine_voxels_percent_of_target_intersections`

100 * five_to_nine_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `ten_or_more_voxels_streamlines`

Target-intersection streamlines visiting at least 10 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `ten_or_more_voxels_percent_of_target_intersections`

100 * ten_or_more_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_most_two_voxels_streamlines`

Target-intersection streamlines visiting at most 2 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_most_two_voxels_percent_of_target_intersections`

100 * at_most_two_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_three_voxels_streamlines`

Target-intersection streamlines visiting at least 3 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_three_voxels_percent_of_target_intersections`

100 * at_least_three_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_five_voxels_streamlines`

Target-intersection streamlines visiting at least 5 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_five_voxels_percent_of_target_intersections`

100 * at_least_five_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_ten_voxels_streamlines`

Target-intersection streamlines visiting at least 10 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: ten_or_more_voxels_streamlines.

### `at_least_ten_voxels_percent_of_target_intersections`

100 * at_least_ten_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: ten_or_more_voxels_percent_of_target_intersections.

### `penetration_p25_unique_target_voxels`

p25 of unique target voxels, using the first histogram bin whose cumulative streamline count reaches the specified percentile.

- Units: voxels.
- Population: target-intersection streamlines.
- Denominator: not applicable.
- Missing/zero: NA when no target-intersection streamlines.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `penetration_median_unique_target_voxels`

median of unique target voxels, using the first histogram bin whose cumulative streamline count reaches the specified percentile.

- Units: voxels.
- Population: target-intersection streamlines.
- Denominator: not applicable.
- Missing/zero: NA when no target-intersection streamlines.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `penetration_p75_unique_target_voxels`

p75 of unique target voxels, using the first histogram bin whose cumulative streamline count reaches the specified percentile.

- Units: voxels.
- Population: target-intersection streamlines.
- Denominator: not applicable.
- Missing/zero: NA when no target-intersection streamlines.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `penetration_p95_unique_target_voxels`

p95 of unique target voxels, using the first histogram bin whose cumulative streamline count reaches the specified percentile.

- Units: voxels.
- Population: target-intersection streamlines.
- Denominator: not applicable.
- Missing/zero: NA when no target-intersection streamlines.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

## endpoint_context

Files or returned table: `*_fibrtraq_endpoint-context.csv (Endpoint coverage)`.

### `context`

reference_grid or a user-named trusted anatomical context mask.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `context_type`

reference_grid checks bounds; binary_mask checks membership in supplied anatomy. These do not establish biological validity.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `valid_streamlines`

Streamlines with valid mapped points, before focal-intersection filtering.

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `endpoint_count`

Available endpoint observations: two per mapped streamline with at least two stored points, one for a single-point record.

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `endpoints_inside`

Endpoints inside this grid/mask.

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `endpoints_outside`

Endpoints not inside this grid/mask; endpoint_count - endpoints_inside.

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `endpoints_nonfinite`

Endpoints whose mapped coordinates are nonfinite (grid row only).

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA in anatomical-mask rows because this is a reference-grid diagnostic.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `endpoints_one_voxel_outside`

Finite endpoints outside bounds with maximum axis excess at most one voxel (grid row only).

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA in anatomical-mask rows because this is a reference-grid diagnostic.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `endpoints_more_than_one_voxel_outside`

Finite endpoints outside bounds with maximum axis excess greater than one voxel (grid row only).

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA in anatomical-mask rows because this is a reference-grid diagnostic.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `streamlines_with_any_endpoint_inside`

Valid mapped streamlines with at least one endpoint inside this grid/mask.

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `streamlines_with_both_endpoints_inside`

Mapped streamlines with all available endpoints inside this grid/mask; a single-point record has one available endpoint.

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `mask_voxels`

Grid volume in voxels or positive anatomical-mask voxels.

- Units: count.
- Population: all valid mapped streamlines/endpoints, not just focal-positive ones.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `endpoints_inside_percent`

100 * endpoints_inside / endpoint_count.

- Units: percent.
- Population: all valid mapped endpoints.
- Denominator: endpoint_count.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `streamlines_with_both_endpoints_inside_percent`

100 * streamlines_with_both_endpoints_inside / valid_streamlines.

- Units: percent.
- Population: all valid mapped streamlines.
- Denominator: valid_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `x_below_grid`

Endpoints below the x grid boundary; distance category: grid . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `x_one_voxel_below`

Endpoints below the x grid boundary; distance category: one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `x_more_than_one_voxel_below`

Endpoints below the x grid boundary; distance category: more than one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `x_above_grid`

Endpoints above the x grid boundary; distance category: grid . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `x_one_voxel_above`

Endpoints above the x grid boundary; distance category: one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `x_more_than_one_voxel_above`

Endpoints above the x grid boundary; distance category: more than one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `y_below_grid`

Endpoints below the y grid boundary; distance category: grid . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `y_one_voxel_below`

Endpoints below the y grid boundary; distance category: one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `y_more_than_one_voxel_below`

Endpoints below the y grid boundary; distance category: more than one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `y_above_grid`

Endpoints above the y grid boundary; distance category: grid . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `y_one_voxel_above`

Endpoints above the y grid boundary; distance category: one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `y_more_than_one_voxel_above`

Endpoints above the y grid boundary; distance category: more than one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `z_below_grid`

Endpoints below the z grid boundary; distance category: grid . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `z_one_voxel_below`

Endpoints below the z grid boundary; distance category: one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `z_more_than_one_voxel_below`

Endpoints below the z grid boundary; distance category: more than one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `z_above_grid`

Endpoints above the z grid boundary; distance category: grid . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `z_one_voxel_above`

Endpoints above the z grid boundary; distance category: one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `z_more_than_one_voxel_above`

Endpoints above the z grid boundary; distance category: more than one voxel . Axis counts overlap when more than one boundary is crossed.

- Units: endpoint count.
- Population: all valid mapped endpoints.
- Denominator: not applicable.
- Missing/zero: 0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `source_path`

Source context mask path; NA for the reference-grid row.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `checksum_md5`

MD5 checksum of the source context mask; NA for the reference-grid row.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `dimensions`

Grid dimensions in x by y by z array order.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

### `affine_match_reference`

Whether the context grid matches the focal/reference affine; the reference-grid row is TRUE.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: all valid mapped streamlines before focal filtering; first discretization candidate.
- Alias of: none.

## local_search_assignments

Files or returned table: `*_fibrtraq_local-search-assignments.csv`.

### `streamline_index`

One-based record position in the input tractogram, not a rank or biological subject identifier.

- Units: streamline index.
- Population: one focal-positive streamline.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

### `start_status`

focal, a target label, no_node_within_radius, ambiguous_nearest_node or endpoint_unavailable. Stored point order defines start/end.

- Units: text.
- Population: one focal-positive streamline endpoint.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

### `start_distance_mm`

Nearest-node scanner-space Euclidean distance from the precise endpoint to a labelled voxel centre; ambiguity can retain a finite distance.

- Units: mm.
- Population: one focal-positive streamline endpoint.
- Denominator: not applicable.
- Missing/zero: 0 means endpoint coincides with a labelled voxel centre; NA when no eligible node/distance is available.
- Method: local_search.
- Alias of: none.

### `start_ambiguous`

Whether competing distinct labels tied or overlapped at the nearest candidate.

- Units: logical.
- Population: one focal-positive streamline endpoint.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

### `end_status`

focal, a target label, no_node_within_radius, ambiguous_nearest_node or endpoint_unavailable. Stored point order defines start/end.

- Units: text.
- Population: one focal-positive streamline endpoint.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

### `end_distance_mm`

Nearest-node scanner-space Euclidean distance from the precise endpoint to a labelled voxel centre; ambiguity can retain a finite distance.

- Units: mm.
- Population: one focal-positive streamline endpoint.
- Denominator: not applicable.
- Missing/zero: 0 means endpoint coincides with a labelled voxel centre; NA when no eligible node/distance is available.
- Method: local_search.
- Alias of: none.

### `end_ambiguous`

Whether competing distinct labels tied or overlapped at the nearest candidate.

- Units: logical.
- Population: one focal-positive streamline endpoint.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

### `paired_target`

Assigned target label only when the two endpoints form an unambiguous focal/target pair.

- Units: text.
- Population: one focal-positive streamline.
- Denominator: not applicable.
- Missing/zero: NA means no assigned pair; consult endpoint statuses.
- Method: local_search.
- Alias of: none.

### `paired_focal_distance_mm`

Distance for the focal-assigned end of a successful pair.

- Units: mm.
- Population: one focal-positive streamline.
- Denominator: not applicable.
- Missing/zero: NA if unpaired; 0 is an exact centre match.
- Method: local_search.
- Alias of: none.

### `paired_target_distance_mm`

Distance for the target-assigned end of a successful pair.

- Units: mm.
- Population: one focal-positive streamline.
- Denominator: not applicable.
- Missing/zero: NA if unpaired; 0 is an exact centre match.
- Method: local_search.
- Alias of: none.

### `paired_max_distance_mm`

Maximum of the two successful-pair endpoint distances.

- Units: mm.
- Population: one focal-positive streamline.
- Denominator: not applicable.
- Missing/zero: NA if unpaired; 0 is an exact centre match.
- Method: local_search.
- Alias of: none.

### `paired_target_intersects_full_path`

Whether the selected voxelized path also occupies an exact voxel of the assigned target. FALSE for an assigned pair means proximity only.

- Units: logical.
- Population: one focal-positive streamline.
- Denominator: not applicable.
- Missing/zero: NA if no pair is assigned.
- Method: local_search.
- Alias of: none.

### `search_radius_mm`

Inclusive configured endpoint search radius, repeated for self-contained rows.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

## local_search_summary

Files or returned table: `tables/*_local-search-summary.csv`.

### `target`

Target receiving at least one successful local-search pair; targets with no pairs have no row.

- Units: text.
- Population: one assigned target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

### `unique_connections`

Successful focal/target pairs; one row per focal-positive streamline in the source assignment table.

- Units: count.
- Population: one assigned target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: local_search.
- Alias of: none.

### `connections_with_full_path_intersection`

Assigned pairs whose selected full path enters the exact assigned target mask.

- Units: count.
- Population: one assigned target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: local_search.
- Alias of: none.

### `connections_without_full_path_intersection`

Assigned pairs whose selected full path does not enter the assigned target mask (proximity only).

- Units: count.
- Population: one assigned target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: local_search.
- Alias of: none.

### `full_path_intersection_percent_of_connections`

100 * connections_with_full_path_intersection / unique_connections.

- Units: percent.
- Population: one assigned target.
- Denominator: unique_connections.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: local_search.
- Alias of: none.

### `focal_distance_median_mm`

median paired_focal_distance_mm; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `target_distance_median_mm`

median paired_target_distance_mm; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `paired_max_distance_q25_mm`

25th percentile paired_max_distance_mm; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `paired_max_distance_median_mm`

median paired_max_distance_mm; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `paired_max_distance_q75_mm`

75th percentile paired_max_distance_mm; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `paired_max_distance_q95_mm`

95th percentile paired_max_distance_mm; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `paired_max_distance_max_mm`

maximum paired_max_distance_mm; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `full_path_target_distance_median_mm`

median paired_target_distance_mm among pairs with full-path intersection; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `full_path_max_distance_median_mm`

median paired_max_distance_mm among pairs with full-path intersection; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `full_path_max_distance_q95_mm`

95th percentile paired_max_distance_mm among pairs with full-path intersection; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `proximity_only_target_distance_median_mm`

median paired_target_distance_mm among proximity-only pairs; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `proximity_only_max_distance_median_mm`

median paired_max_distance_mm among proximity-only pairs; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `proximity_only_max_distance_q95_mm`

95th percentile paired_max_distance_mm among proximity-only pairs; finite distances only, R quantile type 7.

- Units: mm.
- Population: successful assigned pairs in the stated target/subgroup.
- Denominator: not applicable.
- Missing/zero: NA for an empty subgroup; 0 means exact centre match.
- Method: local_search.
- Alias of: none.

### `search_radius_mm`

Inclusive configured endpoint search radius, repeated for self-contained rows.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

## local_search_endpoint_status

Files or returned table: `tables/*_local-search-endpoint-status.csv`.

### `endpoint`

start or end in stored streamline point order.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

### `status`

Endpoint status label as defined in local_search_assignments.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search.
- Alias of: none.

### `streamline_count`

Focal-positive streamlines with this endpoint status.

- Units: count.
- Population: one endpoint/status category.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: local_search.
- Alias of: none.

### `denominator_focal_intersecting_streamlines`

All rows in the focal-positive local-search assignment table.

- Units: count.
- Population: one endpoint/status category.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: local_search.
- Alias of: none.

### `percent_of_focal_intersecting_streamlines`

100 * streamline_count / denominator_focal_intersecting_streamlines.

- Units: percent.
- Population: one endpoint/status category.
- Denominator: denominator_focal_intersecting_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: local_search.
- Alias of: none.

## exploratory_per_voxel

Files or returned table: `ftq_exploratory_summary()$data; explicit research export: *_fibrtraq_exploratory-per-voxel.csv`.

### `fNode_Coords`

Coordinate key x_y_z, using one-based focal-image array indices.

- Units: voxel index key.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Every focal voxel has a key, including zero-assignment voxels.
- Method: selected assignment mode.
- Alias of: none.

### `tracks`

Streamlines allocated to this focal voxel; each streamline contributes at most once to a voxel.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `local`

Allocated streamlines whose entire selected voxelized path is inside the configured local ROI. Local status can overlap target-positive status.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target`

Allocated streamlines that are nonlocal and receive no target under the selected assignment rule.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `x`

One-based focal-image array index along x; not a scanner-space distance.

- Units: voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; zero and NA are invalid.
- Method: focal image grid.
- Alias of: none.

### `y`

One-based focal-image array index along y; not a scanner-space distance.

- Units: voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; zero and NA are invalid.
- Method: focal image grid.
- Alias of: none.

### `z`

One-based focal-image array index along z; not a scanner-space distance.

- Units: voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; zero and NA are invalid.
- Method: focal image grid.
- Alias of: none.

### `leave_with_target`

tracks - local - leave_no_target.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `total_target_hits`

Sum of target associations across this voxel's targets; can exceed tracks.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `target_diversity`

Number of targets with a positive count at this focal voxel.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `primary_target_hits`

Largest target association count at this focal voxel.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `primary_target`

Label with the largest count; ties choose the first target column. This is an association-based label, not proof of anatomical destination.

- Units: text.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: none when no target associations.
- Method: selected assignment mode.
- Alias of: none.

### `primary_target_index`

One-based category code into the primary-target LUT; 0 for no target associations. Do not interpret as a continuous measure.

- Units: category index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no target; positive integers refer to LUT index.
- Method: selected assignment mode.
- Alias of: none.

### `primary_target_fraction`

primary_target_hits / total_target_hits.

- Units: fraction.
- Population: one focal voxel.
- Denominator: total_target_hits.
- Missing/zero: 0 means no numerator events when evaluable; Stored 0 for no target associations, an unevaluable fraction.
- Method: selected assignment mode.
- Alias of: none.

### `target_hits_per_assigned_track`

total_target_hits / tracks.

- Units: fraction.
- Population: one focal voxel.
- Denominator: tracks.
- Missing/zero: 0 means no numerator events when evaluable; Stored 0 for zero tracks; may exceed 1 with multiple targets.
- Method: selected assignment mode.
- Alias of: none.

### `target_entropy_bits`

-sum(p * log2(p)) over positive target proportions p = target_count / total_target_hits.

- Units: bits.
- Population: one focal voxel.
- Denominator: total_target_hits.
- Missing/zero: 0 for a single positive target and also stored 0 for no target associations; use coverage to distinguish them.
- Method: selected assignment mode.
- Alias of: none.

### `target_entropy_fraction`

target_entropy_bits / log2(target_diversity) when diversity > 1.

- Units: fraction.
- Population: one focal voxel.
- Denominator: log2(number of positive targets).
- Missing/zero: Stored 0 for diversity 0 or 1; absence and concentration are distinct.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target_fraction`

leave_no_target / tracks.

- Units: fraction.
- Population: one focal voxel.
- Denominator: tracks.
- Missing/zero: 0 means no numerator events when evaluable; Stored 0 when tracks is zero; unevaluable, not evidence of absence.
- Method: selected assignment mode.
- Alias of: none.

### `leave_with_target_fraction`

leave_with_target / tracks.

- Units: fraction.
- Population: one focal voxel.
- Denominator: tracks.
- Missing/zero: 0 means no numerator events when evaluable; Stored 0 when tracks is zero; unevaluable, not evidence of absence.
- Method: selected assignment mode.
- Alias of: none.

## exploratory_target_totals

Files or returned table: `ftq_exploratory_summary()$target_totals; explicit research export: *_fibrtraq_exploratory-target-totals.csv`.

### `target`

Target label; sorted by decreasing association total, then target label.

- Units: text.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxel_streamline_target_associations`

Sum of this target's per-voxel associations.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxels_with_association`

Focal voxels with a positive count for this target.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `percent_of_focal_voxel_streamline_associations`

100 * focal_voxel_streamline_target_associations / sum of per-voxel tracks.

- Units: percent.
- Population: one target.
- Denominator: sum of per-voxel tracks.
- Missing/zero: 0 means no numerator events when evaluable; Stored 0 for zero total tracks; percentages across targets need not sum to 100.
- Method: selected assignment mode.
- Alias of: none.

## exploratory_lut

Files or returned table: `*_fibrtraq_exploratory-primary-target-lut.csv`.

### `index`

One-based primary_target_index category code; zero is reserved for no target and has no LUT row.

- Units: category index.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `target`

Original target label associated with the category code.

- Units: text.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `color`

Display colour for the target category; a visual key, not a measurement.

- Units: colour code.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

## exploratory_provenance

Files or returned table: `ftq_exploratory_summary()$provenance; explicit research export: *_fibrtraq_exploratory-provenance.csv`.

### `target_intersection_mode`

Primary assignment: full_streamline, endpoint_voxels or local_search. These change target assignment, not the separate full-path endpoint audit.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `local_search_radius_mm`

Inclusive maximum endpoint-to-labelled-voxel-centre distance.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search only.
- Alias of: none.

## run_index

Files or returned table: `fibrtraq-run-index.csv`.

### `subject`

Explicit user-supplied subject identifier; retain leading zeros and suffixes.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_region`

Focal-region label for this run.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `session`

Optional session identifier.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `description`

User-supplied run description.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `prefix`

Filename prefix preserving run identity.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `run_status`

Legacy aggregate state: completed only if the analysis and all requested stages completed without error; failed otherwise. Inspect analysis_status and per-stage fields to identify preserved results.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `biological_qc_status`

Separate anatomical/endpoint QC interpretation; computational success is not biological validation.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `error`

Recorded analysis or presentation error, or empty text if none; inspect individual stage errors when available.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `started_at`

Recorded start timestamp with timezone offset.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `completed_at`

Recorded completion timestamp with timezone offset.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `validation_passed`

Whether executable computational validation passed for the completed run.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `processed_streamlines`

Streamlines processed during the tractogram pass, before focal-positive filtering.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_focal_intersecting_streamlines`

Unique streamlines with at least one selected focal voxel.

- Units: count.
- Population: one run; unique counts are within that run.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `target_intersection_mode`

Primary assignment: full_streamline, endpoint_voxels or local_search. These change target assignment, not the separate full-path endpoint audit.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `selected_track_space`

Coordinate mapping actually used for tractogram points.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `selected_axis_flip`

Comma-separated logical axis flips in x,y,z order.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `local_search_radius_mm`

Inclusive maximum endpoint-to-labelled-voxel-centre distance.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search only.
- Alias of: none.

### `primary_report`

Saved path to main Results HTML ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `exploratory_report`

Saved path to supporting Target report HTML ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `interactive_3d_qc`

Saved path to interactive sampled 3-D QC HTML ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `result_csv`

Saved path to original complete per-voxel table ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `rsc_csv`

Saved path to legacy full RSC table ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `fate_csv`

Saved path to voxel-weighted fate table ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `roi_summary_csv`

Saved path to wide unique-streamline run record ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `endpoint_summary_csv`

Saved path to additive full-path endpoint audit ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `endpoint_penetration_csv`

Saved path to target-penetration histogram ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `endpoint_context_csv`

Saved path to endpoint coverage table ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `local_search_assignments_csv`

Saved path to one-row-per-focal-streamline local-search audit ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `metadata_json`

Saved path to machine-readable run provenance ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `validation_json`

Saved path to machine-readable validation record ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `analysis_status`

analysis stage: completed, failed, not_requested, not_run or not_generated, as applicable. Completed analysis outputs remain available after a presentation failure.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `analysis_error`

Error recorded for the analysis stage; empty text if absent.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `validation_status`

Computational validation stage: passed, failed or not_run. Independent of report-generation status.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `validation_error`

Error recorded for the validation stage; empty text if absent.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `primary_report_status`

primary_report stage: completed, failed, not_requested, not_run or not_generated, as applicable. Completed analysis outputs remain available after a presentation failure.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `primary_report_error`

Error recorded for the primary_report stage; empty text if absent.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `exploratory_report_status`

exploratory_report stage: completed, failed, not_requested, not_run or not_generated, as applicable. Completed analysis outputs remain available after a presentation failure.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `exploratory_report_error`

Error recorded for the exploratory_report stage; empty text if absent.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `interactive_3d_qc_status`

interactive_3d_qc stage: completed, failed, not_requested, not_run or not_generated, as applicable. Completed analysis outputs remain available after a presentation failure.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `interactive_3d_qc_error`

Error recorded for the interactive_3d_qc stage; empty text if absent.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_assignment_mode`

all_focal_voxels counts each occupied focal voxel; first_intersection allocates a streamline only to its first selected focal voxel.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `package_version`

Loaded package identity: version . remote_sha identifies source when present; a development version alone does not.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `package_remote_type`

Loaded package identity: remote type . remote_sha identifies source when present; a development version alone does not.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `package_remote_url`

Loaded package identity: remote url . remote_sha identifies source when present; a development version alone does not.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `package_remote_ref`

Loaded package identity: remote ref . remote_sha identifies source when present; a development version alone does not.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `package_remote_sha`

Loaded package identity: remote sha . remote_sha identifies source when present; a development version alone does not.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `package_built`

Loaded package identity: built . remote_sha identifies source when present; a development version alone does not.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `package_loaded_path`

Loaded package identity: loaded path . remote_sha identifies source when present; a development version alone does not.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

## prepared_inputs

Files or returned table: `fibrtraq-prepared-inputs.csv`.

### `subject`

Explicit user-supplied subject identifier; retain leading zeros and suffixes.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `tractogram`

Path to original tractogram ; unavailable optional inputs do not imply failure.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `reference_source`

Path to source reference image ; unavailable optional inputs do not imply failure.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `bvals_source`

Path to source b-value file (NA when numeric b-values were supplied in memory) ; unavailable optional inputs do not imply failure.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `reference_image`

Path to used/prepared three-dimensional reference ; unavailable optional inputs do not imply failure.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `background_source`

Path to source background image ; unavailable optional inputs do not imply failure.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `background_image`

Path to used/prepared background ; unavailable optional inputs do not imply failure.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `local_mask`

Path to used local-region mask ; unavailable optional inputs do not imply failure.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `reference_source_md5`

MD5 checksum of the corresponding source file, when supplied.

- Units: checksum.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `bvals_source_md5`

MD5 checksum of the corresponding source file, when supplied.

- Units: checksum.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `background_source_md5`

MD5 checksum of the corresponding source file, when supplied.

- Units: checksum.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `reheadered_to_region_grid`

Whether explicit, user-authorized grid-header replacement was requested for prepared images; no interpolation is implied.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

## prepared_volumes

Files or returned table: `prepared/*_prepared-reference-volumes.csv or *_prepared-background-volumes.csv`.

### `source_path`

Original image/b-value source path; bvals_source is NA for numeric b-values supplied in memory.

- Units: file path.
- Population: one selected source volume.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `source_md5`

MD5 checksum of the corresponding source file; NA when no file was supplied.

- Units: checksum.
- Population: one selected source volume.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `bvals_source`

Original image/b-value source path; bvals_source is NA for numeric b-values supplied in memory.

- Units: file path.
- Population: one selected source volume.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `bvals_source_md5`

MD5 checksum of the corresponding source file; NA when no file was supplied.

- Units: checksum.
- Population: one selected source volume.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `volume_index`

One-based source fourth-dimension volume index included in the prepared image.

- Units: volume index.
- Population: one selected source volume.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `b_value`

Diffusion b-value of the selected source volume.

- Units: s/mm^2.
- Population: one selected source volume.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `low_b_max`

Inclusive maximum b-value for low-b averaging.

- Units: s/mm^2.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `source_grid_matches_region_masks`

Whether the source grid already matched regional masks.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `reheadered_to_region_grid`

Whether a prepared derivative received the trusted regional-grid header.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `preparation`

Preparation operation: mean_low_b or the recorded three-dimensional image preparation route.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

## union_inputs

Files or returned table: `prepared/*_regional-union-inputs.csv`.

### `label`

Regional mask label.

- Units: text.
- Population: one union input mask.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `path`

Original regional mask path.

- Units: text.
- Population: one union input mask.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `md5`

MD5 checksum of that source mask.

- Units: text.
- Population: one union input mask.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `nonzero_voxels`

Foreground voxels in this source mask.

- Units: count.
- Population: one union input mask.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

## union_overlaps

Files or returned table: `prepared/*_regional-union-overlaps.csv`.

### `mask_a`

Label of one member of the mask pair.

- Units: text.
- Population: one mask pair.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `mask_b`

Label of one member of the mask pair.

- Units: text.
- Population: one mask pair.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `overlap_voxels`

Foreground voxels shared by both aligned masks; pairs are not independent samples.

- Units: count.
- Population: one mask pair.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

## tractogram_context

Files or returned table: `*_metadata.json: files.context.summary; explicit research export: *_fibrtraq_tractogram-context-summary.csv`.

### `processed_streamlines`

All streamlines visited in the context pass, before focal filtering.

- Units: count.
- Population: whole tractogram context pass.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: context density discretization; independent of focal/target assignment.
- Alias of: none.

### `declared_streamlines`

Tractogram header total; NA if unknown.

- Units: count.
- Population: whole tractogram context pass.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: context density discretization; independent of focal/target assignment.
- Alias of: none.

### `nonzero_tractogram_occupied_voxels`

Reference-grid voxels with at least one sampled or interpolated streamline occupancy.

- Units: count.
- Population: whole tractogram context pass.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: context density discretization; independent of focal/target assignment.
- Alias of: none.

### `total_streamline_voxel_occupancy_events`

Sum of density; one streamline contributes at most once per voxel.

- Units: count.
- Population: whole tractogram context pass.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: context density discretization; independent of focal/target assignment.
- Alias of: none.

### `max_density`

Largest streamline occupancy count at any voxel.

- Units: count.
- Population: whole tractogram context pass.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: context density discretization; independent of focal/target assignment.
- Alias of: none.

### `total_sample_points`

Number of stored points successfully mapped.

- Units: count.
- Population: whole tractogram context pass.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: context density discretization; independent of focal/target assignment.
- Alias of: none.

### `total_voxelization_points`

Number of points after the selected voxelization, including interpolated points when enabled.

- Units: count.
- Population: whole tractogram context pass.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: context density discretization; independent of focal/target assignment.
- Alias of: none.

### `out_of_bounds_sample_points`

Out-of-bounds voxelization points; despite the legacy name this includes interpolated points in segment_interpolation mode.

- Units: count.
- Population: whole tractogram context pass.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: context density discretization; independent of focal/target assignment.
- Alias of: none.

### `mapping_error_streamlines`

Streamlines whose coordinates failed conversion in the context pass.

- Units: count.
- Population: whole tractogram context pass.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: context density discretization; independent of focal/target assignment.
- Alias of: none.

### `tractogram_occupancy_volume_mm3`

Occupied voxel count multiplied by physical voxel volume. This is tractogram occupancy, not anatomical tissue volume.

- Units: mm^3.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `note`

Interpretation reminder that density context is not anatomical brain volume.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_x_min`

min x bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_x_max`

max x bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_y_min`

min y bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_y_max`

max y bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_z_min`

min z bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_z_max`

max z bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_x_min`

min x bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_x_max`

max x bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_y_min`

min y bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_y_max`

max y bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_z_min`

min z bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_z_max`

max z bound of occupied voxels; world bounds use transformed corners of the voxel bounding box.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA when no in-grid occupied voxel.
- Method: selected assignment mode.
- Alias of: none.

## mask_grid

Files or returned table: `ftq_mask_grid_report(..., output_csv = ...) destination`.

### `mask`

Source mask/template/reference image path for this grid audit.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `template_image`

Source mask/template/reference image path for this grid audit.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `reference_mask`

Source mask/template/reference image path for this grid audit.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `mask_dims`

Image array dimensions in x by y by z order.

- Units: voxel dimensions.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `template_dims`

Image array dimensions in x by y by z order.

- Units: voxel dimensions.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `mask_voxel_size`

Raw NIfTI pixdim values in x,y,z order; interpret using header units.

- Units: header spatial units.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `template_voxel_size`

Raw NIfTI pixdim values in x,y,z order; interpret using header units.

- Units: header spatial units.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `dims_match`

Equality of mask and template dimensions.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_size_match`

Equality of raw pixdim values to tolerance 1e-6.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `affine_match`

Equality of unit-normalized declared affines.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `grid_match`

Both dimensions and affine match.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `reference_grid_match_template`

Reference mask shares the template grid.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `nonzero_voxels`

Finite nonzero voxels in the audited mask.

- Units: count.
- Population: one mask/reference comparison.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `overlap_voxels`

Mask/reference foreground intersection, computed only on matching grids.

- Units: count.
- Population: one mask/reference comparison.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `union_voxels`

Mask/reference foreground union, computed only on matching grids.

- Units: count.
- Population: one mask/reference comparison.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `reference_label`

Supplied reference-mask label; NA if no reference was supplied.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `jaccard`

overlap_voxels / union_voxels.

- Units: fraction.
- Population: one mask/reference comparison.
- Denominator: union_voxels.
- Missing/zero: 0 means no numerator events when evaluable; NA if absent references, grids differ, or union is empty.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_x_min`

min x bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_x_max`

max x bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_y_min`

min y bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_y_max`

max y bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_z_min`

min z bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `voxel_bbox_z_max`

max z bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: one-based voxel index.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_x_min`

min x bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: NIfTI world units.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_x_max`

max x bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: NIfTI world units.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_y_min`

min y bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: NIfTI world units.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_y_max`

max y bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: NIfTI world units.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_z_min`

min z bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: NIfTI world units.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

### `world_bbox_z_max`

max z bound of finite nonzero mask voxels, transformed to world coordinates where indicated.

- Units: NIfTI world units.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA for an empty mask.
- Method: selected assignment mode.
- Alias of: none.

## validation_suite

Files or returned table: `ftq_run_validation_suite(..., output_csv = ...) destination`.

### `subject`

Explicit user-supplied subject identifier; retain leading zeros and suffixes.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `session`

Optional session identifier.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `description`

User-supplied run description.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `prefix`

Filename prefix preserving run identity.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `run_status`

Analysis completion state; inspect error and separate report status when present.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `biological_qc_status`

Separate anatomical/endpoint QC interpretation; computational success is not biological validation.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `error`

Recorded analysis or presentation error, or empty text if none; inspect individual stage errors when available.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `started_at`

Recorded start timestamp with timezone offset.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `completed_at`

Recorded completion timestamp with timezone offset.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `validation_passed`

Whether executable computational validation passed for the completed run.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `target_intersection_mode`

Primary assignment: full_streamline, endpoint_voxels or local_search. These change target assignment, not the separate full-path endpoint audit.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `local_search_radius_mm`

Inclusive maximum endpoint-to-labelled-voxel-centre distance.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: local_search only.
- Alias of: none.

### `result_csv`

Saved path to original complete per-voxel table ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `rsc_csv`

Saved path to legacy full RSC table ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `fate_csv`

Saved path to voxel-weighted fate table ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `roi_summary_csv`

Saved path to wide unique-streamline run record ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `endpoint_summary_csv`

Saved path to additive full-path endpoint audit ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `endpoint_penetration_csv`

Saved path to target-penetration histogram ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `endpoint_context_csv`

Saved path to endpoint coverage table ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `metadata_json`

Saved path to machine-readable run provenance ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `validation_json`

Saved path to machine-readable validation record ; paths document the generating workstation; use portable HTML navigation after moving folders.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means not generated or unavailable; consult analysis/report status.
- Method: selected assignment mode.
- Alias of: none.

### `manifest_path`

Manifest, dataset, output or report location recorded for the validation run.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `dataset_root`

Manifest, dataset, output or report location recorded for the validation run.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `output_dir`

Manifest, dataset, output or report location recorded for the validation run.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `report_html`

Manifest, dataset, output or report location recorded for the validation run.

- Units: file path.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `preflight_passed`

Whether input preflight passed; unavailable if not reached.

- Units: logical.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `preflight_errors`

Failed preflight checks of severity error.

- Units: count.
- Population: one validation manifest.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `preflight_warnings`

Failed preflight checks of severity warning.

- Units: count.
- Population: one validation manifest.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `preflight_recommendations`

Failed preflight checks of severity recommendation.

- Units: count.
- Population: one validation manifest.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxels`

Complete focal-voxel count.

- Units: count.
- Population: one validation manifest.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxel_streamline_associations`

Sum of per-voxel tracks.

- Units: count.
- Population: one validation manifest.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `local_associations`

Sum of per-voxel local.

- Units: count.
- Population: one validation manifest.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target_associations`

Sum of per-voxel leave_no_target.

- Units: count.
- Population: one validation manifest.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_with_target_associations`

Sum of tracks - local - leave_no_target.

- Units: count.
- Population: one validation manifest.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `validation_failed_checks`

Semicolon-separated failed computational check names; empty when none.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `endpoint_analysis`

Whether exact_and_dilated1 or none was selected for the additive audit.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxel_streamline_target_associations`

Semicolon-separated target=count association totals; this field is text rather than one numeric total.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

## stat_roi

Files or returned table: `ftq_stat_subject_summary()$roi; saved by the installed study recipe`.

### `subject`

User-supplied subject identifier; preserve leading zeros.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `group`

Optional user-supplied group/condition label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `roi`

Optional focal-region label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxels`

All focal voxels.

- Units: count.
- Population: one subject/focal region; voxel-weighted associations.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `connected_focal_voxels`

Focal voxels with tracks > 0, not just target-positive voxels.

- Units: count.
- Population: one subject/focal region; voxel-weighted associations.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxel_streamline_associations`

Sum of per-voxel tracks.

- Units: count.
- Population: one subject/focal region; voxel-weighted associations.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `local_associations`

Sum of per-voxel local.

- Units: count.
- Population: one subject/focal region; voxel-weighted associations.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target_associations`

Sum of per-voxel leave_no_target.

- Units: count.
- Population: one subject/focal region; voxel-weighted associations.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_with_target_associations`

Sum of tracks - local - leave_no_target.

- Units: count.
- Population: one subject/focal region; voxel-weighted associations.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxel_streamline_target_associations`

Sum of all target columns over all focal voxels.

- Units: count.
- Population: one subject/focal region; voxel-weighted associations.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target_fraction`

leave_no_target_associations / focal_voxel_streamline_associations.

- Units: fraction.
- Population: one subject/focal region.
- Denominator: focal_voxel_streamline_associations.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: selected assignment mode.
- Alias of: none.

### `target_hits_per_focal_voxel_streamline_association`

focal_voxel_streamline_target_associations / focal_voxel_streamline_associations.

- Units: fraction.
- Population: one subject/focal region.
- Denominator: focal_voxel_streamline_associations.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: selected assignment mode.
- Alias of: none.

### `median_target_diversity`

Median number of positive targets over focal voxels with tracks > 0 (includes target-negative assigned voxels).

- Units: targets.
- Population: positive-assignment focal voxels.
- Denominator: not applicable.
- Missing/zero: NA if no assigned voxel; zero means assigned voxels lack targets.
- Method: selected assignment mode.
- Alias of: none.

### `median_entropy_fraction`

Median normalized target entropy over focal voxels with tracks > 0; stored zero for no positive target or only one target is included.

- Units: fraction.
- Population: positive-assignment focal voxels.
- Denominator: not applicable.
- Missing/zero: NA if no assigned voxel; zero can mean no targets or concentrated association.
- Method: selected assignment mode.
- Alias of: none.

## stat_targets

Files or returned table: `ftq_stat_subject_summary()$targets; saved by the installed study recipe`.

### `subject`

User-supplied subject identifier; preserve leading zeros.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `group`

Optional user-supplied group/condition label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `roi`

Optional focal-region label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `target`

Target label.

- Units: text.
- Population: one subject/focal/target row.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxel_streamline_target_associations`

Sum of this target's per-voxel counts.

- Units: count.
- Population: one subject/focal/target row.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `connected_focal_voxels`

Focal voxels with a positive count for this target.

- Units: count.
- Population: one subject/focal/target row.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `target_fraction_of_focal_voxel_streamline_associations`

focal_voxel_streamline_target_associations / sum of per-voxel tracks.

- Units: fraction.
- Population: one subject/focal/target row.
- Denominator: sum of per-voxel tracks.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: selected assignment mode.
- Alias of: none.

### `target_fraction_of_target_hits`

focal_voxel_streamline_target_associations / sum of all target associations.

- Units: fraction.
- Population: one subject/focal/target row.
- Denominator: sum of all target associations.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: selected assignment mode.
- Alias of: none.

## stat_voxels

Files or returned table: `ftq_stat_subject_summary()$voxels; saved by the installed study recipe`.

### `subject`

User-supplied subject identifier; preserve leading zeros.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `group`

Optional user-supplied group/condition label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `roi`

Optional focal-region label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `fNode_Coords`

Coordinate key x_y_z, using one-based focal-image array indices.

- Units: voxel index key.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Every focal voxel has a key, including zero-assignment voxels.
- Method: selected assignment mode.
- Alias of: none.

### `tracks`

Streamlines allocated to this focal voxel; each streamline contributes at most once to a voxel.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `x`

One-based focal-image array index along x; not a scanner-space distance.

- Units: voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; zero and NA are invalid.
- Method: focal image grid.
- Alias of: none.

### `y`

One-based focal-image array index along y; not a scanner-space distance.

- Units: voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; zero and NA are invalid.
- Method: focal image grid.
- Alias of: none.

### `z`

One-based focal-image array index along z; not a scanner-space distance.

- Units: voxel index.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: Positive integer; zero and NA are invalid.
- Method: focal image grid.
- Alias of: none.

### `target_diversity`

Number of targets with a positive count at this focal voxel.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `leave_no_target_fraction`

leave_no_target / tracks.

- Units: fraction.
- Population: one focal voxel.
- Denominator: tracks.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: selected assignment mode.
- Alias of: none.

### `target_hit_events`

Sum of all target associations at this voxel; an alias of exploratory total_target_hits.

- Units: count.
- Population: one focal voxel.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `entropy_fraction`

Normalized target entropy as in exploratory target_entropy_fraction; stored zero for no target or one target.

- Units: fraction.
- Population: one focal voxel.
- Denominator: log2(number of positive targets).
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

## target_unique

Files or returned table: `ftq_target_summary_tables()$unique_streamlines`.

### `target`

Target label.

- Units: text.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `target_family`

Optional user-supplied target family; defaults to unspecified.

- Units: text.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `target_mask_voxels`

Optional supplied target-mask foreground size; NA when not provided.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `unique_target_positive_focal_intersecting_streamlines`

ROI-level unique focal-positive streamlines assigned to this target; read from the source ROI summary, not inferred from voxel totals.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `percent_of_unique_focal_intersecting_streamlines`

100 * unique_target_positive_focal_intersecting_streamlines / unique_focal_intersecting_streamlines from the ROI summary.

- Units: percent.
- Population: one target.
- Denominator: unique_focal_intersecting_streamlines from the ROI summary.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: selected assignment mode.
- Alias of: none.

### `checksum`

Optional user-supplied target-mask checksum.

- Units: checksum.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

## target_associations

Files or returned table: `ftq_target_summary_tables()$voxel_associations`.

### `target`

Target label.

- Units: text.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `target_family`

Optional user-supplied target family; defaults to unspecified.

- Units: text.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `target_mask_voxels`

Optional supplied target-mask foreground size; NA when not provided.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxel_streamline_target_associations`

Sum of this target's per-voxel counts.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `focal_voxels_with_at_least_one_association`

Focal voxels with a positive count for this target.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `percent_of_all_focal_voxel_streamline_associations`

100 * focal_voxel_streamline_target_associations / sum of per-voxel tracks.

- Units: percent.
- Population: one target.
- Denominator: sum of per-voxel tracks.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: selected assignment mode.
- Alias of: none.

## spatial_gradients

Files or returned table: `ftq_stat_spatial_gradients() return table; saved by the installed study recipe`.

### `target`

Target label.

- Units: text.
- Population: one target/axis.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `axis`

Image array axis x, y or z, not a scanner-space anatomical direction.

- Units: text.
- Population: one target/axis.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: selected assignment mode.
- Alias of: none.

### `slope`

Coefficient of the selected array axis in target_count/tracks ~ x + y + z, fitted over finite fractions at positive-assignment focal voxels.

- Units: fraction per voxel index.
- Population: one target/axis.
- Denominator: not applicable.
- Missing/zero: NA when the coefficient is not estimable; unqualified targets have no rows.
- Method: selected assignment mode.
- Alias of: none.

### `p_value`

Nominal linear-model coefficient p-value. Within-subject voxels are not independent biological replicates; this is not a group-level inference.

- Units: probability.
- Population: one target/axis.
- Denominator: not applicable.
- Missing/zero: NA for a non-estimable coefficient.
- Method: selected assignment mode.
- Alias of: none.

### `q_value`

Benjamini-Hochberg adjustment of p_value across returned target-by-axis rows.

- Units: adjusted probability.
- Population: all returned target/axis rows.
- Denominator: not applicable.
- Missing/zero: NA when p_value is unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `connected_voxels`

Positive-assignment focal voxels with a finite target fraction included in the model.

- Units: count.
- Population: one target model.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

### `target_connected_voxels`

Included voxels with a positive target count.

- Units: count.
- Population: one target model.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: selected assignment mode.
- Alias of: none.

## stat_endpoints

Files or returned table: `ftq_stat_subject_summary()$endpoints (optional)`.

### `target`

Target whose exact full-path intersection population is audited.

- Units: text.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `tolerance`

exact or dilated1; the latter includes the 26 neighbouring voxels and is not a uniform millimetre radius.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `tolerance_definition`

Human-readable definition of the endpoint-membership tolerance.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `unique_focal_intersecting_streamlines`

All unique focal-positive streamlines in the run.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_streamlines`

Unique focal-positive streamlines whose selected full path occupies at least one exact target voxel.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_streamlines`

Focal-positive streamlines with at least one focal endpoint at this tolerance, irrespective of target intersection.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_endpoint_streamlines`

Focal-positive streamlines with opposite focal and target endpoints, including ambiguous shared-mask membership and pairs without an exact target intersection.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_with_paired_endpoints`

Exact target-intersection streamlines with opposite focal and target endpoints; this broad count can include shared-mask ambiguity.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_with_any_target_endpoint`

Exact target-intersection streamlines with at least one target endpoint, including ambiguous cases.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_with_any_focal_endpoint`

Exact target-intersection streamlines with any focal endpoint: shared-mask-ambiguous, paired and focal-only classes combined.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_with_target_intersection_streamlines`

Focal-positive streamlines with a focal endpoint and an exact target intersection, including ambiguous membership.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `shared_mask_ambiguous`

Target-intersection streamlines with an endpoint in both focal and target masks at this tolerance; precedence class in the five-class partition.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_focal_target_endpoints`

Unambiguous opposite focal/target endpoints among exact target-intersection streamlines; excludes shared-mask ambiguity.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_endpoint_only`

Unambiguous target endpoint present, no focal endpoint, among target-intersection streamlines.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_only`

Unambiguous focal endpoint present, no target endpoint, among target-intersection streamlines.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `neither_relevant_endpoint`

Neither relevant endpoint in the five-class target-intersection partition.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_mask_voxels`

Focal-mask voxel count at the displayed exact/dilated tolerance.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_mask_voxels`

Target-mask voxel count at the displayed exact/dilated tolerance.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_target_overlap_voxels`

Voxel overlap of focal and target masks at this tolerance.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_percent_of_unique_focal_intersecting_streamlines`

100 * target_intersection_streamlines / unique_focal_intersecting_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: unique_focal_intersecting_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_percent_of_target_intersections`

100 * target_intersection_with_paired_endpoints / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_percent_of_focal_endpoint_streamlines`

100 * paired_endpoint_streamlines / focal_endpoint_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: focal_endpoint_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_endpoint_percent_of_target_intersections`

100 * target_intersection_with_any_target_endpoint / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_percent_of_target_intersections`

100 * target_intersection_with_any_focal_endpoint / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_percent_of_focal_endpoint_streamlines`

100 * focal_endpoint_with_target_intersection_streamlines / focal_endpoint_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: focal_endpoint_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_target_overlap_percent_of_smaller_mask`

100 * focal_target_overlap_voxels / min(focal_mask_voxels, target_mask_voxels).

- Units: percent.
- Population: one target/tolerance.
- Denominator: min(focal_mask_voxels, target_mask_voxels).
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `shared_mask_ambiguous_percent_of_target_intersections`

100 * shared_mask_ambiguous / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_endpoint_only_percent_of_target_intersections`

100 * target_endpoint_only / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_only_percent_of_target_intersections`

100 * focal_endpoint_only / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `neither_relevant_endpoint_percent_of_target_intersections`

100 * neither_relevant_endpoint / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `voxel_size_x_mm`

Reference voxel spacing along x after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `voxel_size_y_mm`

Reference voxel spacing along y after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `voxel_size_z_mm`

Reference voxel spacing along z after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `streamline_voxelization_mode`

sample_points uses stored tractogram points; segment_interpolation densifies consecutive segments before voxelization.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_mode`

Always full_streamline in this additive audit, even when the run's primary assignment uses endpoints or local search.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `subject`

User-supplied subject identifier; preserve leading zeros.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `group`

Optional user-supplied group/condition label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `roi`

Optional focal-region label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

## stat_penetration

Files or returned table: `ftq_stat_subject_summary()$penetration (optional)`.

### `target`

Target whose exact full-path intersection population is audited.

- Units: text.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_streamlines`

Unique focal-positive streamlines whose selected full path occupies at least one exact target voxel.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `one_voxel_touch_streamlines`

Target-intersection streamlines visiting exactly 1 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `one_voxel_touch_percent_of_target_intersections`

100 * one_voxel_touch_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `two_voxel_only_streamlines`

Target-intersection streamlines visiting exactly 2 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `two_voxel_only_percent_of_target_intersections`

100 * two_voxel_only_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `three_to_four_voxels_streamlines`

Target-intersection streamlines visiting 3 or 4 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `three_to_four_voxels_percent_of_target_intersections`

100 * three_to_four_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `five_to_nine_voxels_streamlines`

Target-intersection streamlines visiting 5 through 9 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `five_to_nine_voxels_percent_of_target_intersections`

100 * five_to_nine_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `ten_or_more_voxels_streamlines`

Target-intersection streamlines visiting at least 10 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `ten_or_more_voxels_percent_of_target_intersections`

100 * ten_or_more_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_most_two_voxels_streamlines`

Target-intersection streamlines visiting at most 2 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_most_two_voxels_percent_of_target_intersections`

100 * at_most_two_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_three_voxels_streamlines`

Target-intersection streamlines visiting at least 3 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_three_voxels_percent_of_target_intersections`

100 * at_least_three_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_five_voxels_streamlines`

Target-intersection streamlines visiting at least 5 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_five_voxels_percent_of_target_intersections`

100 * at_least_five_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_ten_voxels_streamlines`

Target-intersection streamlines visiting at least 10 distinct exact target voxels.

- Units: count.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `at_least_ten_voxels_percent_of_target_intersections`

100 * at_least_ten_voxels_streamlines / target_intersection_streamlines.

- Units: percent.
- Population: one target.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `penetration_p25_unique_target_voxels`

p25 of unique target voxels, using the first histogram bin whose cumulative streamline count reaches the specified percentile.

- Units: voxels.
- Population: target-intersection streamlines.
- Denominator: not applicable.
- Missing/zero: NA when no target-intersection streamlines.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `penetration_median_unique_target_voxels`

median of unique target voxels, using the first histogram bin whose cumulative streamline count reaches the specified percentile.

- Units: voxels.
- Population: target-intersection streamlines.
- Denominator: not applicable.
- Missing/zero: NA when no target-intersection streamlines.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `penetration_p75_unique_target_voxels`

p75 of unique target voxels, using the first histogram bin whose cumulative streamline count reaches the specified percentile.

- Units: voxels.
- Population: target-intersection streamlines.
- Denominator: not applicable.
- Missing/zero: NA when no target-intersection streamlines.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `penetration_p95_unique_target_voxels`

p95 of unique target voxels, using the first histogram bin whose cumulative streamline count reaches the specified percentile.

- Units: voxels.
- Population: target-intersection streamlines.
- Denominator: not applicable.
- Missing/zero: NA when no target-intersection streamlines.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `subject`

User-supplied subject identifier; preserve leading zeros.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `group`

Optional user-supplied group/condition label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `roi`

Optional focal-region label.

- Units: text.
- Population: one subject/focal region.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

## target_endpoints

Files or returned table: `ftq_target_summary_tables()$endpoint_streamlines (optional)`.

### `target`

Target whose exact full-path intersection population is audited.

- Units: text.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `tolerance`

exact or dilated1; the latter includes the 26 neighbouring voxels and is not a uniform millimetre radius.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `tolerance_definition`

Human-readable definition of the endpoint-membership tolerance.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `unique_focal_intersecting_streamlines`

All unique focal-positive streamlines in the run.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_streamlines`

Unique focal-positive streamlines whose selected full path occupies at least one exact target voxel.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_streamlines`

Focal-positive streamlines with at least one focal endpoint at this tolerance, irrespective of target intersection.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_endpoint_streamlines`

Focal-positive streamlines with opposite focal and target endpoints, including ambiguous shared-mask membership and pairs without an exact target intersection.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_with_paired_endpoints`

Exact target-intersection streamlines with opposite focal and target endpoints; this broad count can include shared-mask ambiguity.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_with_any_target_endpoint`

Exact target-intersection streamlines with at least one target endpoint, including ambiguous cases.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_with_any_focal_endpoint`

Exact target-intersection streamlines with any focal endpoint: shared-mask-ambiguous, paired and focal-only classes combined.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_with_target_intersection_streamlines`

Focal-positive streamlines with a focal endpoint and an exact target intersection, including ambiguous membership.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `shared_mask_ambiguous`

Target-intersection streamlines with an endpoint in both focal and target masks at this tolerance; precedence class in the five-class partition.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_focal_target_endpoints`

Unambiguous opposite focal/target endpoints among exact target-intersection streamlines; excludes shared-mask ambiguity.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_endpoint_only`

Unambiguous target endpoint present, no focal endpoint, among target-intersection streamlines.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_only`

Unambiguous focal endpoint present, no target endpoint, among target-intersection streamlines.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `neither_relevant_endpoint`

Neither relevant endpoint in the five-class target-intersection partition.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_mask_voxels`

Focal-mask voxel count at the displayed exact/dilated tolerance.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_mask_voxels`

Target-mask voxel count at the displayed exact/dilated tolerance.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_target_overlap_voxels`

Voxel overlap of focal and target masks at this tolerance.

- Units: count.
- Population: one target/tolerance.
- Denominator: not applicable.
- Missing/zero: 0 means no counted events; NA means unavailable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_percent_of_unique_focal_intersecting_streamlines`

100 * target_intersection_streamlines / unique_focal_intersecting_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: unique_focal_intersecting_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_percent_of_target_intersections`

100 * target_intersection_with_paired_endpoints / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `paired_percent_of_focal_endpoint_streamlines`

100 * paired_endpoint_streamlines / focal_endpoint_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: focal_endpoint_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_endpoint_percent_of_target_intersections`

100 * target_intersection_with_any_target_endpoint / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_percent_of_target_intersections`

100 * target_intersection_with_any_focal_endpoint / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_percent_of_focal_endpoint_streamlines`

100 * focal_endpoint_with_target_intersection_streamlines / focal_endpoint_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: focal_endpoint_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_target_overlap_percent_of_smaller_mask`

100 * focal_target_overlap_voxels / min(focal_mask_voxels, target_mask_voxels).

- Units: percent.
- Population: one target/tolerance.
- Denominator: min(focal_mask_voxels, target_mask_voxels).
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `shared_mask_ambiguous_percent_of_target_intersections`

100 * shared_mask_ambiguous / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_endpoint_only_percent_of_target_intersections`

100 * target_endpoint_only / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `focal_endpoint_only_percent_of_target_intersections`

100 * focal_endpoint_only / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `neither_relevant_endpoint_percent_of_target_intersections`

100 * neither_relevant_endpoint / target_intersection_streamlines.

- Units: percent.
- Population: one target/tolerance.
- Denominator: target_intersection_streamlines.
- Missing/zero: 0 means no numerator events when evaluable; NA when denominator is zero.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `voxel_size_x_mm`

Reference voxel spacing along x after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `voxel_size_y_mm`

Reference voxel spacing along y after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `voxel_size_z_mm`

Reference voxel spacing along z after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute.

- Units: mm.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `streamline_voxelization_mode`

sample_points uses stored tractogram points; segment_interpolation densifies consecutive segments before voxelization.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_intersection_mode`

Always full_streamline in this additive audit, even when the run's primary assignment uses endpoints or local search.

- Units: text.
- Population: one run.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `target_family`

Optional user-supplied target family; defaults to unspecified.

- Units: text.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.

### `checksum`

Optional user-supplied target-mask checksum.

- Units: checksum.
- Population: one target.
- Denominator: not applicable.
- Missing/zero: NA means unavailable or not applicable.
- Method: additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius.
- Alias of: none.
