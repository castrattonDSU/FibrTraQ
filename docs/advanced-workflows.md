# Advanced workflows and installation

[Back to the quick start](../README.md) · [Interpretation guide](interpretation.md) · [S1c mouse-study recipe](mouse-study-workflow.md)

Start with the bundled example and one-analysis template in the README. This
page retains the configuration, provenance and lower-level workflows needed
for a reproducible study or an integration with other software.

**On this page:** [Fixed versions](#install-a-fixed-version) ·
[Private SSH access](#private-repository-access-with-ssh) ·
[Verify an installed version](#verify-an-installed-version) ·
[Input checklist](#data-checklist) · [Regional study](#run-a-regional-study) · [Manifest](#manifest-workflow) ·
[Core function](#core-counting-function) · [Grid safeguards](#coordinate-and-grid-safeguards) ·
[Outputs](#what-fibrtraq-writes) · [Older report links](#sharing-and-opening-reports) ·
[Function reference](#main-functions) · [Performance](#backend-and-performance)

## Install a fixed version

The README installs the current version on `main`. For a study,
choose a reviewed release tag or full 40-character commit SHA, then use the
same value on every workstation. Replace the placeholder below with that
actual identifier; a package version alone does not identify a development
snapshot.

```r
install.packages("remotes")
approved_ref <- "<RELEASE_TAG_OR_40_CHARACTER_COMMIT_SHA>"
remotes::install_github(
  "castrattonDSU/FibrTraQ",
  ref = approved_ref,
  dependencies = NA,
  upgrade = "never",
  force = TRUE
)
```

Restart R after installing. For private access, accept the repository invitation
and use your configured GitHub credential. Never store an authentication token
in this script or an output report.

## Private repository access with SSH

If SSH access is already configured, install the same approved revision with:

```r
approved_ref <- "<APPROVED_40_CHARACTER_COMMIT_SHA>"
remotes::install_git(
  "git@github.com:castrattonDSU/FibrTraQ.git",
  ref = approved_ref,
  git = "external",
  dependencies = NA,
  upgrade = "never",
  force = TRUE
)
```

GitHub's [SSH connection guide](https://docs.github.com/en/authentication/connecting-to-github-with-ssh)
explains the one-time setup. Restart R after installation. For an independent
review, use the installed package; `devtools::load_all()` instead loads a source
checkout.

## Verify an installed version

Verify the loaded package and the main public interfaces before running study
data:

```r
library(FibrTraQ)

required_exports <- c(
  "ftq_run_files",
  "ftq_run_regional_study",
  "ftq_run_example",
  "ftq_local_search_tables"
)

required_run_files_arguments <- c(
  "assignment_mode",
  "search_radius_mm"
)

stopifnot(
  all(required_exports %in% getNamespaceExports("FibrTraQ")),
  all(required_run_files_arguments %in%
        names(formals(FibrTraQ::ftq_run_files)))
)

description <- packageDescription("FibrTraQ")
cat(
  "FibrTraQ version: ", description$Version, "\n",
  "Loaded from: ", find.package("FibrTraQ"), "\n",
  "GitHub SHA: ", if (is.null(description$RemoteSha)) {
    "not recorded"
  } else {
    description$RemoteSha
  }, "\n",
  sep = ""
)
```

For an exact-commit review, compare the printed GitHub SHA with the approved
SHA. Retain that identity with the run record; a development version number by
itself does not identify an exact source state.

## Data checklist

### Required files

| File | Requirement |
|---|---|
| Tractogram | One TrackVis `.trk` file containing 3-D streamline points |
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

The `.trk` coordinate convention must also be compatible with that image grid.
If the convention is unknown, the package can compare plausible mappings for
a pilot run, but the selected mapping and 3-D QC still require review.

## Run a regional study

Use `ftq_run_regional_study()` when each subject has the same named set of
regional masks and each region should take a turn as the focal ROI. For every
subject and requested focal region, all other regions become targets.

For example, 6 subjects with 4 focal regions produce 24
subject-by-focal-region analyses. They remain 6 independent biological
subjects, not 24 independent samples.

The function:

- validates subject names, mask labels, non-empty masks, and image grids, and
  reports nonbinary values for review;
- optionally prepares a 3-D mean low-b reference from a 4-D image;
- builds a documented same-grid union of the regional masks by default;
- runs each focal region against all remaining regions;
- writes standard and Target reports, maps, tables, and validation files;
- creates one interactive 3-D viewer per run; and
- creates one study-level HTML index.

### Regional-study input structure

| Argument | R structure |
|---|---|
| `tractograms` | Named vector or list with one `.trk` path per subject |
| `region_masks` | Named outer list; each subject contains the same named regional paths |
| `reference_images` | Optional named vector/list with one image per subject |
| `bvals` | Named paths or numeric vectors, required for 4-D references |
| `endpoint_contexts` | Optional named outer list of trusted context masks |
| `assignment_mode` | One primary target-assignment mode for every run |
| `search_radius_mm` | Positive physical radius used by `"local_search"` |
| `output_dir` | One study-level output directory |

Subject names must match across every supplied object. Regional labels must
match across subjects, although different subjects may have different image
dimensions.

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

study <- ftq_run_regional_study(
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
assessment. A technically valid run can still require biological review.
Anatomical endpoint containment is unavailable unless trusted context masks
are supplied.

### Regional local-mask behavior

The default `local_masks = "union"` builds a same-grid union of all supplied
regions for each subject. Under that setting, `local` means that a streamline
stays inside the regional union; it does not mean that the streamline stays
inside the focal region.

| Setting | Meaning |
|---|---|
| `local_masks = "union"` | Default: build one same-grid union per subject |
| `local_masks = NULL` | Use the active focal mask as the local mask |
| Named subject paths | Use independently validated local masks |
| `focal_regions = c("cerebellum", "pontine")` | Run only selected focal regions |

For 4-D references, supply a named `bvals` path or numeric vector for every
subject. Use a distinct output directory when testing alternative inputs.

Open `fibrtraq-report-index.html` first. It links every standard report,
Target report, and interactive viewer. The `exploratory_report` index field
and `write_exploratory_report` argument retain their existing names for
compatibility. The two CSV indexes preserve the
run status, output paths, input checksums, prepared references, backgrounds,
and local masks.

## Manifest workflow

Use a YAML or JSON manifest when an analysis must be shared, reviewed, or
rerun without editing R code. Paths may be relative to `dataset_root`.

Minimal YAML example:

```yaml
dataset_root: /absolute/path/to/dataset
subject: subject1
description: cerebellum

tractogram: tractography/subject1.trk
fnode: rois/subject1_cerebellum.nii.gz
ref_image: anatomy/subject1_reference.nii.gz
local:

targets:
  cerebrum: rois/subject1_cerebrum.nii.gz
  pontine: rois/subject1_pontine.nii.gz
  medullaspinal: rois/subject1_medullaspinal.nii.gz

analysis:
  target_intersection_mode: full_streamline
  local_search_radius_mm: 2

outputs:
  output_dir: derivatives/fibrtraq-v0.3.0
  overwrite: false

report:
  background_image: anatomy/subject1_reference.nii.gz
```

`ref_image` defines coordinate conversion between the tractogram and voxel
grid. `report.background_image` supplies visual anatomy. They serve different
roles even when they point to the same file.

Preflight and run the manifest:

```r
library(FibrTraQ)

preflight <- ftq_preflight("/absolute/path/to/manifest.yml")
stopifnot(preflight$passed)

run <- ftq_run(
  "/absolute/path/to/manifest.yml",
  backend = "auto"
)
stopifnot(run$validation$passed)

browseURL(run$report$html)

viewer <- ftq_interactive_qc_viewer(
  run,
  overwrite = FALSE,
  open = TRUE
)
```

The installed example manifest is also available at:

```r
system.file(
  "examples",
  "manifest-example.yml",
  package = "FibrTraQ"
)
```

`ftq_preflight()` checks input paths, target labels, coordinate and target
assignment options, mask values, image dimensions and affines, output
collisions, optional background compatibility, strict endpoint-context
requirements, and expected regression baselines before the tractogram pass.

`ftq_run()` writes the standard report and core files. Use the viewer call
above for interactive QC. Use `ftq_run_files()` or
`ftq_run_regional_study()` when the standard report, Target report,
interactive viewer, and report index should all be generated automatically.

## Core counting function

`FibrTraQ()` is the exported scientific counting engine used by every standard
workflow. The higher-level functions do not implement alternative counting
methods. They prepare inputs, call the same core function, validate its result,
and organize the output files:

```text
ftq_run_files()
  -> ftq_run_regional_study()
    -> ftq_run()
      -> fibrtraq()
        -> FibrTraQ()
```

The lowercase `fibrtraq()` function is a compatibility alias that forwards its
arguments unchanged to `FibrTraQ()`. It contains no separate analysis logic.

Call `FibrTraQ()` directly when the ROI volumes and tractography are already
loaded in R or when a program needs direct control over the in-memory
calculation:

```r
library(FibrTraQ)

focal <- freesurferformats::read.nifti1.data(
  "/absolute/path/to/focal.nii.gz"
)
local <- freesurferformats::read.nifti1.data(
  "/absolute/path/to/local.nii.gz"
)
targets <- list(
  cerebrum = freesurferformats::read.nifti1.data(
    "/absolute/path/to/cerebrum.nii.gz"
  ),
  pontine = freesurferformats::read.nifti1.data(
    "/absolute/path/to/pontine.nii.gz"
  ),
  medullaspinal = freesurferformats::read.nifti1.data(
    "/absolute/path/to/medullaspinal.nii.gz"
  )
)
tracks <- freesurferformats::read.dti.trk(
  "/absolute/path/to/tractogram.trk"
)
reference_image <- "/absolute/path/to/reference.nii.gz"

result <- FibrTraQ(
  fNode = focal,
  localROI = local,
  cNode_list = targets,
  cNode_labels = names(targets),
  all_tracks = tracks,
  track_space = "auto",
  ref_image = reference_image,
  backend = "auto",
  assignment_mode = "full_streamline",
  search_radius_mm = 2
)

per_voxel <- result
roi_summary <- attr(result, "fibrtraq_roi_summary")
rsc <- rsc_summary(
  result,
  percent = TRUE,
  include_legacy_tracks = FALSE
)
endpoint <- ftq_endpoint_tables(result)
```

`result` contains one row per focal voxel and retains the ROI-level,
coordinate-mapping, and endpoint results as attributes. The principal
in-memory outputs are:

| Object | Meaning |
|---|---|
| `result` | Per-focal-voxel streamline and target-association table |
| `attr(result, "fibrtraq_roi_summary")` | ROI-level unique-streamline counts and run totals |
| `rsc_summary(result, ...)` | Focal-voxel/streamline/target association and RSC summary |
| `ftq_endpoint_tables(result)` | Endpoint-pair, target-penetration, context, and endpoint-metadata tables |

This call uses the current scientific defaults: `"any"` discretization,
all occupied focal voxels, stored streamline sample points, full-path target
assignment, and exact plus one-voxel endpoint analysis. The explicit mode and
radius make the analysis contract visible; the radius is used only by
`"local_search"`. Set these arguments explicitly when running a sensitivity
analysis or reproducing an older method.

A direct core call does not write CSV files, NIfTI maps, validation JSON,
reports, or the interactive 3-D viewer. Use `ftq_run_files()`,
`ftq_run_regional_study()`, or `ftq_run()` when those reproducible workflow
artifacts are required. Direct callers are also responsible for file-level
preflight and for proving the tractogram's coordinate convention. `FibrTraQ()`
rejects local and target arrays with different dimensions, rejects different
unit-normalized NIfTI affines when both objects declare a valid qform/sform,
and applies the same checks to an explicit reference image whenever that image
can affect mapping or local search. Coordinate-matrix ROIs are checked for finite,
integer-compatible, in-bounds voxel indices, but they have no affine metadata;
their spatial provenance remains the caller's responsibility.

The automatic coordinate search shown above is an initial diagnostic, not an
anatomical alignment test. When `axis_flip` is omitted, the search still scores
candidate mappings, but stops before counting if the best candidate reflects
any axis. Sampled ROI overlap cannot distinguish a correct orientation from a
reflection that happens to intersect more of the focal mask. Independently
verify the tractogram convention and its alignment with the masks, then set
`track_space` and `axis_flip` explicitly. Validated reflections remain supported;
no flips is not a universal solution. Review the selected mapping and 3-D QC
before interpreting results, and pin `discretization` for a production analysis.
See `?FibrTraQ` for the complete core interface.

## Coordinate and grid safeguards

These checks matter for every workflow:

1. **Use one subject grid at a time.** Focal, local, target, reference,
   background, and endpoint-context images must have compatible dimensions and
   affines. File workflows check these before counting, and the direct
   `FibrTraQ()` engine independently rejects incompatible local, target, and
   used reference grids.
2. **Treat automatic tract mapping as a pilot.** `track_space = "auto"` and
   `discretization = "any"` compare plausible mappings. If `axis_flip` is omitted
   and the best candidate reflects an axis, the analysis stops for an explicit
   orientation decision. An unreflected candidate also needs independent
   alignment review; a successful run does not establish anatomical validity.
3. **Set known coordinates explicitly.** Supported track spaces include
   `"voxel"`, `"voxelmm"`, `"vox10"`, `"world"`, and `"trackvis"`. Once verified,
   pin both `track_space` and `axis_flip`, including explicit `FALSE` values for
   axes that must not be reflected.
4. **Do not repair alignment by changing headers casually.**
   `reheader_images_to_region_grid = TRUE` is only for data whose voxel-index
   correspondence was independently proven and whose source header is known
   to be non-authoritative.
5. **Use trusted-affine regridding conservatively.** `ftq_regrid_mask()` cannot
   repair a mask whose spatial metadata are wrong.
6. **Review the whole streamline.** A target-positive streamline can extend
   outside the displayed foreground or anatomical mask. The viewer shows the
   sampled tractography points; it does not clip them to anatomy.
7. **Distinguish density context from anatomy.** Without a reference image,
   report heatmaps use a grayscale background derived from the tractogram under
   the selected mapping. This background cannot independently validate that
   mapping. Use trusted anatomical images and spatial provenance to assess
   alignment.

Use `ftq_mask_grid_report()` for a read-only grid audit:

```r
grid <- ftq_mask_grid_report(
  mask = "/absolute/path/to/mask.nii.gz",
  template_image = "/absolute/path/to/focal.nii.gz",
  output_csv = "/absolute/path/to/mask-grid-report.csv"
)

grid
```

When several trusted regional masks already share the focal grid,
`ftq_make_union_mask()` is preferable to resampling an untrusted
whole-brain mask.

## What FibrTraQ writes

Open the report index, then **Results**. See the [main interpretation summary](../README.md#read-the-main-results) for the RSC formula and counting populations. The main report leads with the
RSC summary (RSC last), core maps, complete voxel data and the relevant
fate/coverage interpretation. **Destination** names the target or local
comparison, and **Total number of fNode voxels** states the denominator.
Detailed endpoint QC, target patterns,
provenance and diagnostics remain accessible in supporting sections.

Reports link to a self-contained `OUTPUTS_README.html` and a Markdown
copy defining columns, maps, units, counting populations, denominators, missing
values, aliases and coordinate conventions. The same
[output dictionary is included with the package](../inst/output-guide/OUTPUTS_README.md).
The Target report reuses the shared guide instead of writing another copy.
If an older or customized guide already exists, it is retained and the report
links a separate current guide with a content-specific filename.

The high-level workflows organize the following outputs under the chosen directory:

| Output | Interpretation |
|---|---|
| Report index | First page to open; one Results action per run and separate stage status |
| Standard HTML report | RSC summary, core maps, complete data and fate/coverage, with readable Methods and Validation details |
| Target HTML report | Primary-target, diversity, entropy, and many-target spatial summaries; shares its `figures/` folder with five images |
| Primary-target lookup CSV | Matches integer target codes to labels; saved at the focal-region root and readable in Excel |
| Interactive 3-D QC | Target-balanced quota subset of positive streamlines, ROI extents, end voxels, and optional advanced hit/class controls |
| Original per-voxel CSV | Complete focal-voxel/streamline association counts with the preserved canonical schema |
| All-voxel CSV and column-role CSV | The same complete counts with explicit one-based image-array coordinates; the companion file distinguishes coordinates from targets with colliding names |
| Local-fraction CSV and maps | `local / tracks`; zero-assignment voxels are unevaluable. Use the companion map validity mask |
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

The affected part of each new run's layout is:

```text
subject/
└── focal-region/
    ├── <prefix>_fibrtraq_exploratory-primary-target-lut.csv
    ├── figures/
    │   ├── <prefix>_fibrtraq_target-report.html
    │   └── [five target-report images]
    ├── maps/
    │   ├── [core NIfTI maps]
    │   └── <prefix>_fibrtraq_tractogram-density.nii.gz  [when needed]
    └── reports/
        └── <run>/
            ├── [main results HTML]
            ├── report-figures/
            └── tables/
                └── <prefix>_ranked_<metric>_voxels.csv
```

Core scientific tables and their schemas are unchanged. A subject-level
`prepared/` folder is created only when a workflow must save a derived image
or mask, such as a low-b reference, reheadered reference, or regional union.
New standard runs do not create `context/` or `exploratory-report/` folders.
Existing results are not reorganized or deleted. Compatibility aliases remain
in the full CSVs and are identified in the dictionary.

The standard Target report retains its HTML, five images, and primary-target
lookup table. Extra exploratory per-voxel, target-total and provenance CSVs,
and exploratory NIfTI maps are reserved for explicit research exports.
Use `ftq_exploratory_summary(result_csv, output_dir = research_dir,
template_image = focal_mask, write_maps = TRUE)` with a separate, new research
directory to obtain those exports. Calling the summary with
`output_dir = NULL, write_maps = FALSE` returns the calculations in memory.
These derived metrics retain their method and denominator qualifications;
they are not additional validated measures of biological connectivity.

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

When no anatomical background is supplied, `.trk` workflows can generate an
exact tractogram-density image for report orientation. This is visualization
context only. It is not MRI, a brain mask, tissue volume, or evidence of
anatomical alignment. The density image is stored with the core maps;
its summary and settings are retained in the run-metadata JSON rather than
separate context files.

Generated metadata preserve source paths, checksums, user-supplied subject
labels, the normalized target-assignment mode, and the local-search radius for
reproducibility. Review these fields and use de-identified subject labels
before sharing an output bundle publicly.

## Sharing and opening reports

Copy or ZIP the **whole results folder**, keeping its subfolders together.
Extract the ZIP before opening `fibrtraq-report-index.html` in a browser.
The summary uses relative links, so the complete folder can move between
Windows, macOS and Linux without changing its links. Spaces in folder names
are supported. CSV provenance retains the original input and execution paths.

For an older summary whose Windows links contain `%5C`, install the corrected
package and open the repair script in RStudio:

```r
file.edit(system.file(
  "scripts", "fibrtraq-rebuild-report-index.R", package = "FibrTraQ"
))
```

Set its `report_dir` to the folder containing `fibrtraq-run-index.csv`, then
run the script. Open the new `fibrtraq-report-index-portable.html`. This rebuilds
only the summary, preserves the original files, and requires no tractography
rerun. It also works after moving a complete results folder to another machine.

## Main functions

| Function | Purpose |
|---|---|
| `ftq_run_example()` | Run and verify the bundled known-answer analysis |
| `ftq_run_files()` | Run one complete focal-versus-target analysis from files |
| `ftq_run_regional_study()` | Run all requested subject-by-focal regional analyses |
| `ftq_preflight()` | Validate a manifest before the tractogram pass |
| `ftq_run()` | Run one manifest-defined analysis |
| `FibrTraQ()` | Core in-memory scientific counting engine used by every workflow |
| `ftq_interactive_qc_viewer()` | Build the rotatable 3-D QC viewer |
| `ftq_endpoint_tables()` | Retrieve endpoint, penetration, context, and metadata tables |
| `ftq_target_summary_tables()` | Separate unique-streamline and voxel-association target summaries |
| `ftq_mask_grid_report()` | Audit dimensions, affines, and optional mask overlap without modifying files |

Open a function's help page or the complete package index in R:

```r
?ftq_run_files
?ftq_run_regional_study
?ftq_run
?FibrTraQ
help(package = "FibrTraQ")
```

## Backend and performance

The high-level file runners use `backend = "auto"` by default: they select the
compiled C++ backend when available and otherwise fall back to the reference R
backend. The manifest example passes this setting explicitly. Direct
`ftq_run()` and `FibrTraQ()` calls retain the R backend default for backward
compatibility.

Large `.trk` files are streamed in chunks by default. No GPU, CUDA, Python,
ANTs, or other external analysis program is required to run the R package.
Use `backend = "r"` only for validation or debugging.
