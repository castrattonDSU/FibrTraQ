# Installation, study settings and output files

[README](../README.md) · [Interpretation guide](interpretation.md) · [Spatial verification](spatial-verification.md) · [S1c mouse recipe](mouse-study-workflow.md)

## Install a fixed version

Use the same reviewed tag or full commit for every workstation. Choose a revision containing `fibrtraq_single()` and `fibrtraq_study()`. For a revision already on GitHub, replace the example reference below with its tag or full commit SHA.

```r
install.packages("remotes")
remotes::install_github(
  "castrattonDSU/FibrTraQ", ref = "<REVIEWED_TAG_OR_COMMIT>",
  dependencies = NA, upgrade = "never"
)
```

For a private repository, use an account with access and its existing GitHub credential. Keep authentication tokens out of scripts. To review unpublished changes locally, open the checkout containing `DESCRIPTION` in RStudio and run `remotes::install_local(".", dependencies = NA, upgrade = "never")`. Building from source requires a suitable compiler, such as Rtools on Windows or the Xcode Command Line Tools on macOS.

## Verify an installed version

Restart R and run:

```r
library(FibrTraQ)
stopifnot(setequal(getNamespaceExports("FibrTraQ"),
                   c("fibrtraq_single", "fibrtraq_study")))
packageDescription("FibrTraQ")[c("Version", "RemoteSha")]
```

GitHub installations record their revision in `RemoteSha`. Keep the local source revision separately when installing a local checkout.

## Data checklist

### Required files

| File | Requirement |
|---|---|
| Tractogram | One TrackVis `.trk` or MRtrix `.tck` file containing 3-D streamline points |
| Focal ROI | One non-empty 3-D NIfTI mask (`.nii` or `.nii.gz`); values should be finite and binary |
| Target ROIs | One or more named, non-empty 3-D NIfTI masks; values should be finite and binary |
| Output directory | A new or intentionally replaceable directory for derived files |

### Recommended or optional files

| File | Why it helps |
|---|---|
| 3-D reference image | Coordinate conversion and spatial context for reports/QC |
| 4-D diffusion image plus matching b-values | FibrTraQ can prepare a documented mean low-b reference |
| Anatomical background image | Clearer static reports and 3-D orientation |
| Local-region mask | Defines the region used for the non-exclusive `local` classification |
| Brain, brain-plus-cord, or other context mask | Audits whether streamline endpoints lie inside trusted anatomy |

For each subject, all active masks and images must represent the same voxel
grid. Matching dimensions alone are not enough; NIfTI affines must also be
compatible. FibrTraQ rejects grid mismatches rather than silently resampling
or reheadering them.

Focal, local, and target masks use nonzero values as foreground; preflight
reports nonbinary values for review. Endpoint context masks are stricter and
must be finite, binary, three-dimensional, and on the focal grid.

The tractogram's coordinates must also be compatible with that image grid.
For `.trk`, if the convention is unknown, the package can compare plausible
mappings for a pilot run, but the selected mapping and 3-D QC still require
review. Native `.tck` has a defined world-coordinate convention and follows
the requirements below.

### MRtrix TCK inputs

Use a `.tck` path in the same `tractogram` argument accepted by
`fibrtraq_single()`, or in the subject’s `tractograms` entry for
`fibrtraq_study()`. No conversion to `.trk` is required. The report, target figures, density maps, endpoint outputs and interactive 3-D viewer
use the same analysis as the `.trk` workflow.

| Property | Native TCK support |
|---|---|
| Storage | Uncompressed, single-file `.tck` with an embedded binary payload; `.tck.gz` and external payload files are not supported |
| Numeric types | `Float32LE`, `Float32BE`, `Float64LE` and `Float64BE` |
| Reading | `track_reader = "auto"` streams by default; `"stream"` and `"memory"` are also available |
| Coordinates | RAS+ world coordinates in millimeters, mapped through the NIfTI affine |
| Coordinate arguments | `track_space = "auto"` or `"world"`; `axis_flip = NULL` or `c(FALSE, FALSE, FALSE)` |
| Weights | Ordinary streamline counts; external SIFT2 weights are not applied |

For TCK, `"auto"` selects world coordinates deterministically. It does not
search for the orientation with the most focal-mask hits. Do not add axis
flips or override the tractogram affine/voxel size. A negative axis direction
or oblique orientation encoded in the NIfTI affine is valid and is handled by
that affine.

The mapping reference must be a genuine NIfTI image with declared spatial
units (millimeters, meters or micrometers) and a valid, invertible qform or
sform. When both transforms are declared, they must agree within 0.001 mm at
the voxel-center grid corners. Missing units, absent or invalid transforms,
and conflicting qform/sform geometry stop the run. Resolve those source
metadata with the data producer rather than changing headers to pass a check.
FibrTraQ converts the declared image units to millimeters for coordinate
mapping; it does not register or resample the inputs.

If `reference_image` is omitted, the focal NIfTI supplies the coordinate
reference. A trusted, same-grid anatomical image remains strongly recommended
for independent spatial review.

TCK support allows human research tractograms to use these workflows. It does
not establish anatomical validity for every dataset or change the meaning of
FibrTraQ's unweighted counts and RSC measures. Review each dataset using the
[spatial-verification guide](spatial-verification.md).

## Run a regional study

Use `fibrtraq_study()` when each subject has the same named set of
regional masks and each region should take a turn as the focal ROI. For every
subject and requested focal region, all other regions become targets.

For example, 6 subjects with 4 focal regions produce 24
subject-by-focal-region analyses. They remain 6 independent biological
subjects, not 24 independent samples.

The study validates each subject’s grids and labels. It can prepare a mean low-b reference and a union of regional masks. Its single HTML report includes maps, tables, target patterns and a 3D viewer for each selected analysis.

### Regional-study input structure

Use a CSV sheet for routine cohorts. Copy the [two-subject template](../inst/examples/study-sheet.csv) into the study folder and enter one row per subject. The installed copy is available at `system.file("examples", "study-sheet.csv", package = "FibrTraQ")`.

| Sheet column | What to enter |
|---|---|
| `subject` | A unique text ID, such as `001` |
| `tractogram` | The subject's `.trk` or `.tck` file |
| `reference_image` | Optional same-subject reference. Leave a cell blank for a subject without one |
| `mask:<region>` | One path per subject to that region's NIfTI mask. Provide at least two regional columns |

Relative input paths are resolved from the directory containing the CSV. The `output_dir` argument follows normal R path conventions. Define the tractograms, regional masks, and optional references together in the sheet. To supply them directly in R, use the named-input form below.

| Argument | R structure |
|---|---|
| `tractograms` | Named vector or list with one path per subject |
| `region_masks` | Named outer list with the same regional names for each subject |
| `reference_images` | Optional named vector or list of reference paths |

Both forms accept `focal_regions`, `assignment_mode`, and advanced options. For example, `bvals` supplies named paths or numeric vectors for 4-D references, and `endpoint_contexts` supplies trusted context masks. Each subject's images must share its grid. Different subjects may have different image dimensions.

### Check the inputs before counting

Add `check_only = TRUE` to either public call with the settings you intend to use. The result contains `passed` for the whole request, `index` for subject-by-focal-region status, and `checks` for individual findings. The check leaves the requested output folder untouched. Correct input errors and verify anatomical alignment before running the same analysis with `check_only = FALSE`, the default.

### Complete regional-study example

```r
library(FibrTraQ)

data_root <- "/absolute/path/to/study-data"
results_dir <- "/absolute/path/to/fibrtraq-regional-results"
s1 <- file.path(data_root, "subject1")
s2 <- file.path(data_root, "subject2")

tractograms <- c(
  subject1 = file.path(s1, "subject1.trk"),
  subject2 = file.path(s2, "subject2.trk")
)

region_masks <- list(
  subject1 = c(
    cerebellum = file.path(s1, "subject1_cerebellum.nii.gz"),
    cerebrum = file.path(s1, "subject1_cerebrum.nii.gz"),
    pontine = file.path(s1, "subject1_pontine.nii.gz"),
    medullaspinal = file.path(s1, "subject1_medullaspinal.nii.gz")
  ),
  subject2 = c(
    cerebellum = file.path(s2, "subject2_cerebellum.nii.gz"),
    cerebrum = file.path(s2, "subject2_cerebrum.nii.gz"),
    pontine = file.path(s2, "subject2_pontine.nii.gz"),
    medullaspinal = file.path(s2, "subject2_medullaspinal.nii.gz")
  )
)

reference_images <- c(
  subject1 = file.path(s1, "subject1_reference.nii.gz"),
  subject2 = file.path(s2, "subject2_reference.nii.gz")
)

study <- fibrtraq_study(
  tractograms = tractograms,
  region_masks = region_masks,
  reference_images = reference_images,
  output_dir = results_dir,
  description = "regional-connectivity"
)

stopifnot(
  all(study$index$run_status == "completed"),
  all(study$index$validation_passed)
)

browseURL(study$report_index)
```

With two subjects and four regions, this call creates eight analyses. The
study index records every subject, focal region, selected coordinate mapping,
validation result, and report/output path:

```r
study$index[, c(
  "subject",
  "focal_region",
  "run_status",
  "validation_passed",
  "biological_qc_status",
  "primary_report",
  "exploratory_report",
  "interactive_3d_qc"
)]
```

`validation_passed` is the computational gate.
`biological_qc_status` is a separate `pass`, `review`, or `unavailable`
assessment. A computationally valid run can still require anatomical review.
Anatomical endpoint containment is unavailable unless trusted context masks
are supplied.

### Regional local-mask behavior

Both study and single-subject runs use the active focal mask as the local region by default. This means that a streamline is local if all its selected path voxels stay inside the focal region. To reproduce a regional study run under the 0.3.5 study default, set `local_masks = "union"`. Under that rule, a local streamline can leave the focal region and stay within the union of all supplied regional masks. The choice affects local and fate summaries.

| Setting | Meaning |
|---|---|
| `local_masks = NULL` | Default. Use the active focal mask as the local mask |
| `local_masks = "union"` | Build one same-grid union of the subject’s regions |
| Named subject paths | Use independently validated local masks |
| `focal_regions = c("cerebellum", "pontine")` | Run only selected focal regions |

For 4-D references, supply a named `bvals` path or numeric vector for every
subject. Use a distinct output directory when testing alternative inputs.

Open `study$report_index` and select the subject and focal region. The run
and prepared-input CSVs record status, paths, input checksums and any derived
images. The `write_exploratory_report` and `write_qc_viewer` arguments control
the target and viewer sections of the report.

## Coordinate and grid safeguards

Both workflows check dimensions and NIfTI geometry before counting. Review each subject’s masks and tractogram on an independent anatomical reference, then keep the supported `track_space`, `axis_flip`, discretization and voxelization with the run record. Native TCK uses the NIfTI affine directly. TRK automatic mapping stops if its preferred candidate would require an axis flip. A density background generated from the tractogram can show the selected geometry but cannot independently verify it. [Spatial verification](spatial-verification.md)

## What FibrTraQ writes

The console prints the report path and the folder to share. Pass `open = TRUE` to open the HTML automatically, or reopen it at `study$report_index`. The RSC destination table and heatmaps appear near the top. Each target name links to its NIfTI map when present, and **Total number of fNode voxels** gives the RSC denominator. Detailed checks, methods and the interactive 3D viewer are accessible in the same report. [Interpretation summary](../README.md#read-the-main-results)

Refer to the [package field dictionary](../inst/output-guide/OUTPUTS_README.md)
for field definitions, units and column aliases.

The high-level workflows organize the following outputs under the chosen directory:

| Output | Interpretation |
|---|---|
| HTML report | RSC, spatial previews, fate/coverage, targets, methods, validation and interactive QC for the selected analysis |
| Primary-target lookup CSV | Maps integer target codes to labels in `tables/` |
| Interactive 3-D QC | Target-balanced quota subset of positive streamlines, ROI extents, end voxels, and optional advanced hit/class controls |
| Original per-voxel CSV | Complete focal-voxel/streamline association counts with the preserved canonical schema |
| All-voxel CSV and column-role CSV | The same complete counts with explicit one-based image-array coordinates; the companion file distinguishes coordinates from targets with colliding names |
| Local-fraction CSV and figure | `local / tracks` when tracks is positive. Unassigned focal voxels have an undefined fraction |
| ROI-summary CSV | ROI-level unique-streamline totals and run metadata |
| RSC CSV | Target associations and fraction of focal voxels connected to each target |
| Fate CSV | Local, leave-with-target, and leave-without-target association groups |
| Report connection-summary CSV | Compact local/target rows with unique streamline connections, connected focal voxels, RSC and the total focal-voxel denominator; its existing file schema is preserved |
| Report fate-summary CSV | Mutually exclusive unique-streamline partition into local, leave-with-cNode, and leave-without-cNode categories |
| End-voxel summary CSV | Exact and one-voxel end-voxel membership and paired ends |
| End-voxel penetration CSV | Histogram of unique target voxels visited per target-positive streamline |
| Endpoint coverage CSV (`*_endpoint-context.csv`) | Reference-grid and optional anatomical-mask containment across valid mapped streamlines, before focal filtering |
| Local-search assignment CSV | One row per focal-positive streamline, recording end-voxel labels and distances, the assigned cNode, and whether the voxelized path also enters that cNode |
| NIfTI maps | Voxel-resolved result fields on the focal grid |
| Metadata and validation JSON | Machine-readable methods, input provenance, arithmetic, reconciliation, regression and separate biological-QC status |
| Ranked CSVs (`*_ranked_<metric>_voxels.csv`) | Every focal voxel, including zero counts, ranked by the named target or metric, with deterministic coordinate tie-breaking |

For one analysis the HTML is in the focal-region directory. For several
subject/region combinations a single study report is written at the root.

```text
output/
├── fibrtraq-run-index.csv
├── fibrtraq-prepared-inputs.csv
└── subject/
    └── focal-region/
        ├── figures/    PNGs and 3D viewer data
        ├── tables/     Scientific CSVs, rankings and lookup tables
        ├── maps/       Core NIfTIs and any density background
        └── ...         Analysis metadata and validation JSON
```

A subject's `prepared/` directory contains derived inputs only when needed.

Use the original `*_per-voxel-hits.csv` for general analysis. The all-voxel
export is for inspection and coordinate lookup; do not treat its added
coordinate columns as target counts when loading it with `read.csv()`.

If a report or viewer fails after the analysis finishes, the run index retains
the completed scientific files and records the failing stage. Inspect
`analysis_status`, `validation_status` and the individual report statuses:
the legacy aggregate `run_status` can be `"failed"` even when the analysis
completed. A missing optional input, an unavailable anatomical QC assessment,
and a computation failure are different conditions. Actionable warnings and
failures remain visible.

When no anatomical background is supplied, `.trk` and `.tck` workflows can generate an
exact tractogram-density image for report orientation. This is visualization
context only. It is not MRI, a brain mask, tissue volume, or evidence of
anatomical alignment. The density image is stored with the core maps;
its summary and settings are retained in the run-metadata JSON rather than
separate context files.

Generated metadata preserve source paths, checksums, user-supplied subject
labels, the normalized target-assignment mode, and the local-search radius for
reproducibility. Review these fields and use de-identified subject labels
before sharing an output bundle publicly.

## Windows output paths

Use a short writable output root, for example `C:/FibrTraQ/run01`, and concise
subject, description, focal-region and target labels. Subject, region and table filenames contribute to the full path length.

Before regional analysis, FibrTraQ checks planned Windows output paths. If a
path exceeds the supported length, the run stops with that path and guidance
to choose a shorter destination. Resolve the error before rerunning; existing
output folders are not renamed or deleted. The same analysis settings and
input masks can be reused with a new output directory.

This protects the portable report workflow; it does not establish support for
arbitrarily long paths. Some Windows file operations and browser tools retain
a 260-character full-path limit even when other applications can create the
file. Keep complete results folders together and use a short location when
moving reports between workstations.

## Sharing and opening reports

Send the whole directory labeled **Share this folder**. A single-run HTML sits two levels below that directory. If sending a ZIP, extract it before opening the HTML. Select the subject and focal region in a multi-analysis report. The report’s figures, tables and 3D viewer use relative paths. The metadata contain the original input paths.

## Method reference

Use `?fibrtraq_single` and `?fibrtraq_study` for the complete arguments and return values. The [interpretation guide](interpretation.md) explains target assignment, focal-voxel coverage and endpoint QC.
