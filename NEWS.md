# FibrTraQ changes

## 0.3.4 (CRAN resubmission candidate; unreleased)

- Check planned Windows output paths before regional analysis and explain how
  to shorten an unsupported destination. Existing output names and scientific
  calculations are preserved.
- Shortened generated temporary output names in the bundled example and tests
  to accommodate deeper temporary directories such as CRAN's Windows setup.
- Replaced README links to files excluded from the source archive with public
  repository URLs, preserving the linked sections.
- Added release checks for deeper Windows temporary paths, README file links
  and the additional incoming spelling diagnostics. Unexpected check notes
  require review before submission.

## 0.3.3 (first CRAN submission; not accepted)

- Clarified supported tractogram formats, outputs and spatial-verification
  requirements in the package description for the first CRAN submission.
- Declared optional plotting/development dependencies used by the installed
  public-data validation scripts and the minimum testing-library version.
- Isolated reader and figure-capture test substitutes from other packages'
  namespaces, preserving the scientific regression assertions.

- Added native MRtrix `.tck` inputs to the existing file, manifest and regional
  workflows, including streamed and in-memory reading of single-file
  Float32/Float64 data in either byte order. Standard and Target reports,
  density context, endpoint outputs, maps and the interactive 3-D QC viewer
  use the same scientific counting engine as `.trk` inputs.
- TCK coordinates follow their RAS+ world-millimeter convention through the
  NIfTI affine without orientation guessing. TCK runs require declared image
  units and valid, consistent NIfTI transforms, and reject extra flips or
  incompatible coordinate overrides. Existing `.trk` mapping behavior is
  preserved. Inputs are not automatically registered or resampled, and
  streamline counts remain unweighted.
- Documented TCK input requirements and the distinction between format support
  and independent anatomical validation in the README and workflow guides.
- Preserved focal-image voxel sizes in endpoint audit tables when no separate
  reference image is supplied. Non-isotropic TRK inputs previously displayed
  unit voxel sizes in those metadata fields; streamline counts and mappings
  are unchanged.
- Wrapped long endpoint-figure target labels without truncating identifiers,
  preventing overlap with the explanatory footnote.
- Added pinned public native-TCK validation scripts and an independent reader
  oracle, with a separate audit of Float32 conversion boundary effects.
- Kept in-memory TRK reading compatible with both traditional track lists and
  the compact track collections returned by recent `freesurferformats`
  versions. Native stored coordinates are requested explicitly when supported.

## 0.3.1 (research preview candidate; unreleased)

- Aligned package and citation metadata for the next research preview. The
  proposed release tag is `v0.3.1-rc1`; the candidate must complete the release
  checks before that tag or a public release is created.
- Added a practical spatial-verification guide for checking independent
  anatomical evidence, recording the selected mapping, and reusing verified
  settings. Clarified the single-subject README instructions.
- Updated the release checklist to separate automated package checks, actual
  report/browser review, and the public source/history review.
- Added a native Windows browser acceptance workflow for the installed
  synthetic example, including relocated report links, downloads, 3-D
  interactions, and screenshots for visual review. Its results are recorded
  separately from anatomical validation.
- Corrected initial 3-D scene framing to enclose diagonal paths after rotation,
  including in narrow windows. Display coordinates and scientific results are
  unchanged; user-controlled zoom remains available.

## 0.3.0.9002

- Automatic coordinate selection now stops before counting when `axis_flip`
  is omitted and the best candidate reflects any axis. A sampled ROI-overlap
  score can favor an anatomically incorrect reflection, so users must verify
  the mapping and explicitly set the axis flips. Explicit, independently
  validated reflections remain supported.
- Clarified how to pin a verified mapping and that tractogram-derived grayscale
  backgrounds are density context, not independent anatomical validation.

## 0.3.0.9001

- Improved the 3-D QC viewer with free trackball rotation, edge dragging to
  roll, and two-finger pinch zoom. Reset restores the original anatomical view
  and zoom. Drawing updates now follow the screen refresh cycle and reuse the
  viewport and rotation transform, reducing work during a drag. Sampling,
  scientific results, display layers, and the initial orientation are unchanged.

- Relabeled the main report's connection summary as **RSC summary**, its
  denominator as **Total number of fNode voxels**, and its displayed target
  column as **Destination**, without changing the calculation or CSV schema.
- Simplified newly generated public outputs: the **Target report** and its
  five images now share the focal-region `figures/` folder, and the primary-target
  lookup CSV sits at the focal-region root. Extra exploratory CSVs, NIfTI maps,
  and duplicated output guides are no longer generated by the standard report
  workflow. Explicit research exports through `ftq_exploratory_summary()`
  remain available, and existing output folders are preserved.
- Moved generated tractogram-density context images into the core `maps/`
  folder and consolidated context settings and summaries into run metadata.
  Empty `prepared/` folders are no longer created; necessary prepared images
  and masks are retained.
- Renamed the main report's image folder to `report-figures/`. Ranked CSVs now
  contain every focal voxel, including zero-count rows, and use
  `*_ranked_<metric>_voxels.csv` filenames instead of top-50 filenames.

- Changed the license for the current source and package-authored synthetic
  fixture from MIT to the GNU General Public License, version 2 or later
  (`GPL-2.0-or-later`). Earlier copies, including `v0.3.0-rc2`, retain their
  original MIT terms. Third-party dependency licenses are unchanged.
- Aligned package, citation, contributor and fixture notices. The licensing
  change does not alter analysis code or scientific reference values.

## 0.3.0 (release candidate v0.3.0-rc2)

- Simplified first use to GitHub installation, the bundled example, and the
  Results report. Added a portable mouse-study recipe with exact input filenames
  and separate guidance for advanced workflows and interpretation.
- Labeled RSC as the percentage of focal voxels connected and explained its
  all-focal-voxel denominator and overlapping rows. RSC values and scientific
  output schemas are unchanged; these percentages need not sum to 100%.
- Reject Windows-reserved subject/region folder names and trailing periods
  before regional outputs are created, preserving portability without silently
  changing sample identifiers.
- Fixed the installed DTOR start-to-finish recipe's exploratory-report call
  to use the public argument names and retain the completed run's method and
  radius. A regression exercises that exact call with the bundled fixture.

- Focused the main Results report on the connection summary, core maps,
  complete per-voxel data and explicit fate/coverage populations. Detailed
  exploratory, endpoint, validation and provenance records remain accessible;
  existing scientific CSV schemas, map values and generation defaults are
  retained. RSC is displayed last, and endpoint counts and percentages use
  explicitly matched populations.
- Added complete coordinate-enriched voxel exports with a companion column-role
  dictionary, preserving genuine targets named x, y or z. Added local-fraction
  CSV and NIfTI views with an explicit validity mask for zero-assignment voxels.
  Top-50 tables remain optional ranked previews and no longer repeat the
  leave-no-target column. Existing files are preserved.
- Added self-contained HTML and Markdown output guides covering field
  definitions, units, denominators, missing/zero values, aliases, maps and
  assignment/tolerance dependence. Reusing an identical guide is safe; older
  or customized guides are retained beside a content-specific current guide.
- Regional run indexes now retain completed analysis outputs after later report
  or viewer failures and distinguish analysis, validation and presentation
  stages. Report entry points use explicit subject labels and support optional
  one-based axial preview slice selection.
- Corrected exploratory labels to distinguish voxel associations from unique
  streamlines and the configured local region from the focal region. CSV-based
  primary and exploratory reports recover saved method/radius provenance when available
  and visibly identify unverified fallback metadata. Documented target ties and
  separate assignment/target coverage when reporting exploratory distributions.
- Clarified the example's evidence: seven exact scientific-table comparisons
  and 18 independent scientific claims are distinct from the four nonempty
  HTML-file checks. Browser interaction and anatomical validity remain separate
  validation requirements.

- Fixed report-summary links and browser opening on Windows. Report links use
  relative forward-slash URLs and preserve spaces, reserved characters and
  literal percent signs in filenames when the complete results folder moves.
  Figure creation also preserves literal percent signs in output-folder names.
  Standalone report downloads now resolve from the report's actual directory.
  Added an installed script to rebuild an existing summary from its run index
  and existing HTML files without rerunning the analysis.

- Development builds now use a distinct version and record the loaded package
  path plus available GitHub installation fields in run metadata and regional
  run indexes. This makes source identity auditable without recording access
  tokens or other credentials.

- Added `ftq_run_files()` as the beginner-oriented entry point for one complete
  analysis. Users supply a tractogram path, focal-mask path, named target-mask
  paths, and an output directory; the function writes the standard and
  exploratory reports, endpoint outputs, maps, interactive 3-D QC viewer, and
  a report index without requiring a manifest or manual object loading.
- Added `ftq_run_regional_study()` as the production high-level entry point for
  multi-subject regional analyses. A single call validates named file inputs,
  prepares 3-D low-b references from 4-D images, builds regional union masks,
  runs each focal region against all remaining targets, and writes standard
  reports, exploratory reports, endpoint outputs, 3-D QC viewers, and a
  study-level report index.
- Fixed explicit non-resampled reheadering of 3-D reference images so stale
  source-array header attributes cannot override the regional-mask affine.
  When the mask template omits spatial units and the source declares valid
  units, the prepared derivative now preserves that source unit declaration.
- Added selectable primary target assignment through the direct-call
  `assignment_mode` argument (saved as `target_intersection_mode` in manifests
  and provenance): backward-compatible `"full_streamline"`
  (all selected path voxels), `"endpoint_voxels"`, or `"local_search"`.
  Local search uses precise endpoints, an inclusive
  `search_radius_mm = 2` default, and Euclidean scanner-space distance to
  labelled voxel centres. NIfTI metre and micrometre spatial units are
  converted to millimetres; unknown units are recorded and treated as
  millimetres with a warning. Competing distinct-label ties and overlapping
  labels are left ambiguous. Paper-terminology aliases are normalized in
  provenance.
  These modes operationalize the broad strategies compared by Yeh et al.
  (2019) without claiming exact reproduction of unspecified implementation
  details. Focal qualification/allocation remains in the centralized
  `FibrTraQ()` engine, and the existing exact/dilated endpoint analysis remains
  an additive full-path audit.
- Added a reconciled local-search end-voxel audit without changing established
  result-table schemas. Local-search runs now retain one row per focal-positive
  streamline with endpoint labels, scanner-space distances, the assigned
  cNode, and whether the selected voxelized path also enters that cNode. The
  new `ftq_local_search_tables()` helper, dedicated assignment CSV, validation
  checks, report tables, and distance/reconciliation figure distinguish direct
  target intersections from proximity-only local-search connections.
- Added default additive endpoint-versus-intersection analysis in the existing
  tractogram pass. New result attributes and dedicated CSVs report exact and
  one-voxel-sensitivity focal-target endpoint pairs, five mutually exclusive
  endpoint classes, and exact sparse target-penetration histograms without
  changing existing per-voxel, RSC, fate, ROI-summary, or map schemas.
- Added optional trusted binary endpoint context masks with strict
  dimension/affine preflight checks, reference-grid boundary diagnostics, and a
  separate biological QC status that does not redefine computational success.
- Added endpoint classification, sensitivity, penetration, and context sections
  to an advanced end-voxel QC section in static reports, plus exact/one-voxel
  modes, class filters, and start/end markers in the existing
  selected-target-assignment-positive 3D QC subset.
- Refined the standard report for publication-facing review. It now uses a
  compact preflight and validation summary, adds a local row to the RSC table,
  reports unique-streamline connection and mutually exclusive fate summaries
  from ROI-level counters, explains every primary figure, and moves detailed
  end-voxel and density-distribution diagnostics into clearly labeled advanced
  sections. Existing per-voxel, RSC, fate, ROI-summary, and map schemas remain
  unchanged.
- Added display-only robust scaling to spatial density heatmaps so isolated
  high-density voxels do not flatten visible structure. Local and cNode density
  rows use distinct palettes and independently labeled scales; source CSV and
  NIfTI values are unchanged.
- Reoriented the interactive 3D QC viewer to an anatomical RAS display basis
  with superior up and a rotating R/L, A/P, S/I compass. End voxels are visible
  by default, dense hit markers are initially hidden, and uncommon controls are
  grouped under an advanced disclosure. Viewer text now states that the
  target-quota subset is accepted in tractogram scan order and is neither
  random nor prevalence weighted.
- Added `ftq_endpoint_tables()` and opt-in endpoint additions to subject- and
  target-summary helpers. `endpoint_analysis = "none"` retains a regression and
  performance path with all established metrics unchanged.
- Added `ftq_tractogram_context()` and automatic report/QC fallback to exact
  tractogram-derived density context when no anatomical background image is
  supplied. The context is labeled as visualization-only and is not treated as
  anatomical brain volume.
- Added mask-grid harmonization utilities: `ftq_mask_grid_report()` for
  non-destructive dimension/affine diagnostics, `ftq_make_union_mask()` for
  building aligned biological union masks from same-grid regional ROIs, and
  conservative `ftq_regrid_mask()` for trusted-affine nearest-neighbor mask
  regridding.
- Preflight now distinguishes dimension mismatch from same-dimension
  affine/grid mismatch and recommends union masks when aligned regional masks
  are available.
- Direct and manifest-driven runs now reject local, target, and used reference
  images that do not share the focal voxel grid. NIfTI affines are compared
  after normalizing declared spatial units when both images declare a valid
  qform/sform; unavailable transforms are reported without treating RNifti's
  fallback identity as provenance. Coordinate-matrix ROIs are checked for
  finite, integer-compatible, in-bounds voxel indices without assuming affine
  metadata they do not contain.
- Added a compiled Rcpp counting backend selected by `backend = "auto"` when
  available, with reference R parity tests and profiling support for benchmark
  runs.
- Upgraded HTML report figures with full-slice context plus focal-zoom panels
  for tract counts, target maps, leave-no-target maps, and exploratory maps.
- Upgraded the interactive 3D QC viewer with full focal ROI extent, target ROI
  extent clouds, target-hit markers, segment display modes, sampled streamline
  controls, and tractogram-derived context when anatomy is absent.
- Added `ftq_example_data()` and `ftq_run_example()` with a package-owned,
  checksummed TrackVis/NIfTI known-answer fixture. One call now writes the
  complete reports and 3-D QC output, then requires exact reproduction of seven
  versioned scientific tables and 18 independently declared scientific values.
- Refocused the README on installation, required inputs, the verified example,
  and beginner-to-advanced workflows. Completed evidence, acceptance gates,
  biological-validity limits, and template-only public benchmark candidates
  are now distinguished in `VALIDATION.md` and the extended validation docs.
- Added explicit `focal_assignment_mode` support. The recommended default is now
  `all_focal_voxels`, which counts every distinct focal voxel occupied by a
  discretized stored streamline sample point and is invariant to streamline
  point-order reversal.
- Added explicit `streamline_voxelization_mode = "sample_points"` metadata and
  manifest/API support. Under this mode FibrTraQ does not rasterize every
  continuous line segment between stored points.
- Retained the legacy one-focal-voxel behavior as
  `focal_assignment_mode = "first_intersection"` for compatibility and
  sensitivity analyses.
- Added an ROI-level unique-streamline summary table with focal, local,
  leave-no-target, leave-with-target, local-by-target exclusive fate groups,
  target, association-count, coordinate-mapping, and mode metadata. Unique
  streamline totals should be read from this summary rather than inferred by
  summing per-voxel association counts.
- Added `ftq_target_summary_tables()` to separate ROI-level unique target
  streamlines from focal-voxel/streamline/target association totals and avoid
  denominator ambiguity.
- Harmonized `discretization = "any"` so one selected voxelization governs
  focal, local, and target classification for each streamline.
- Added preflight diagnostics for focal/local, focal/target, local/target, and
  target/target mask overlap. Local is reported as a non-exclusive flag.
- Updated automatic coordinate mapping selection so target/local ROI hits are
  reported for audit but not used to select or tie-break mappings.
