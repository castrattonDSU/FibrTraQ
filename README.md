# FibrTraQ

[![R-CMD-check](https://github.com/castrattonDSU/FibrTraQ/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/castrattonDSU/FibrTraQ/actions/workflows/R-CMD-check.yaml)

FibrTraQ is an analysis package that **assesses and quantifies the connectivity of a specified focal node (fNode) to an unlimited set of connecting nodes (cNodes)**. The primary outputs are in two formats: voxel-wise fNode heatmaps for each defined cNode, and computed region-to-region spatial connectivity (RSC) percentage values in summary tables that describe how much of the fNode is connected to each cNode. The package also generates a 3D quality control viewer. FibrTraQ is suitable to quantify how much of a fNode has connections with each target cNode, inspect the spatial distribution of those connections, and examine how full-path or endpoint-based target assignment affects the result. Of note, interpretation is reliant on the plausibility of tractography, accuracy of node delineations and spatial harmonization across file inputs.

To run the package, you need a TrackVis `.trk` or MRtrix `.tck` tractogram and binary NIfTI masks (`.nii` or `.nii.gz`) for the fNode and each cNode. A same-subject, aligned reference image (such as b0, T1, T2 or FA) is optional for execution and strongly recommended for spatial review.

**Research preview:** FibrTraQ's tests establish the documented counting behavior and computational reproducibility. Each dataset still requires independent spatial verification before its results are interpreted. Follow the [spatial-verification guide](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/spatial-verification.md) when bringing new inputs into the package.

**Start here:** [Install and try the example](#install-and-try-the-rudimentary-example) ·
[Single-subject analysis](#run-fibrtraq-on-a-single-subject) · [Read the results](#read-the-main-results) ·
[Multiple subjects](#run-a-regional-study) · [Troubleshooting](#common-problems)

<a id="1-install-fibrtraq"></a>

## Install and try the rudimentary example

FibrTraQ requires **R 4.0.0 or later**. In a fresh R or RStudio session, install
the current version on GitHub's `main` branch:

```r
install.packages("remotes")
remotes::install_github(
  "castrattonDSU/FibrTraQ",
  ref = "main",
  dependencies = NA,
  upgrade = "never",
  force = TRUE
)
```

Version 0.3.4 is being prepared for CRAN resubmission after the 0.3.3
incoming checks identified Windows report-path and documentation-link issues.
FibrTraQ is not yet available through `install.packages("FibrTraQ")`; use the
GitHub installation above until a CRAN release is announced.

<a id="2-verify-the-installation"></a>

Below is a supplementary basic example for observing the logic of the package. 
To run this, first **Restart R**, then run the rudimentary example. No external data are needed:

```r
library(FibrTraQ)

demo <- ftq_run_example(
  output_dir = file.path(tempdir(), "ftq-example"),
  overwrite = TRUE
)

stopifnot(demo$validation$passed)
browseURL(demo$report_index)
```

On the page that opens, choose **Results**. The main report contains an RSC summary table (which may be difficult to interpret in this basic example), an RSC histogram for visual examination, aggregate fNode streamline density heatmaps, fNode streamline connection density heatmaps per cNodes, maps of fNode voxels with no streamlines, fNode local connection density heatmaps, and a fate summary table.

At the top of the page are links to the **Target report** and a 3-D quality control viewer. Example files are in R's temporary directory; use a persistent output folder for your own analyses.

In the 3-D viewer, drag to rotate freely, or drag near an edge to roll the view.
Scroll or pinch with two fingers to zoom. **Reset Orientation** restores the
initial anatomical view and zoom; the compass follows the rotation.

<details>
<summary>Compilers and reproducible versions</summary>

The compiled backend is built during installation. If R requests a compiler, use Apple Command Line Tools on macOS, Rtools matching your R version on Windows, or an R-compatible C++ toolchain on Linux.

`main` changes over time. For a study, install a reviewed release tag or full commit SHA and use that same revision across workstations. See [fixed-version installation and verification](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#install-a-fixed-version).

</details>

## Run FibrTraQ on a single subject

Start with a tractogram, one fNode and at least one specified cNode. Use nonempty, finite, binary 3-D NIfTI masks (`.nii` or `.nii.gz`). All masks and images must share the subject's voxel grid, including compatible NIfTI affines. The template below includes a 3-D reference image for spatial review. If it is unavailable, omit `reference_image`; the resulting density background cannot independently verify anatomy.

**This is a template for your own files.** The filenames below are illustrative; replace each one with the name of your existing input. Set the four directories, subject identifier and ROI labels for your dataset. A publication describing the package and validation against a prior manual simplified version of this quantitative analysis (Cottam et al., 2023) is in preparation. Contact nickcottam4@gmail.com with questions regarding this dataset.

```r
library(FibrTraQ)

subject_id <- "subject-001"
node_dir <- "/absolute/path/to/input-nodes"
trk_dir <- "/absolute/path/to/input-trk"
ref_dir <- "/absolute/path/to/ref-image" # same subject and spatial grid as the masks
output_dir <- "/absolute/path/to/new-fibrtraq-results"

targets <- c(
  cerebrum = file.path(node_dir, "cerebrum.nii.gz"),
  pontine = file.path(node_dir, "pontine.nii.gz"),
  medullaspinal = file.path(node_dir, "medullaspinal.nii.gz")
)

analysis <- ftq_run_files(
  tractogram = file.path(trk_dir, "tractogram.trk"),
  focal_mask = file.path(node_dir, "cerebellum.nii.gz"),
  target_masks = targets,
  reference_image = file.path(ref_dir, "reference.nii.gz"),
  output_dir = output_dir,
  subject = subject_id,
  focal_label = "cerebellum",
  assignment_mode = "full_streamline"
)

analysis$index[, c(
  "subject", "focal_region", "analysis_status", "validation_status"
)]
browseURL(analysis$report_index)
```

For an MRtrix tractogram, replace `"tractogram.trk"` with your `.tck` filename
in the same call. The reports, maps and 3-D viewer use the same workflow.
TCK points are read as RAS+ world coordinates in millimeters and mapped through
the NIfTI affine; leave `track_space = "auto"` or use `"world"`, with no extra
axis flips. The NIfTI reference (or focal mask when no reference is supplied)
must have valid, consistent spatial geometry and declared units. See the
[TCK input requirements](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#mrtrix-tck-inputs) for details.

On Windows, choose a short output directory such as `C:/FibrTraQ/run01` and
keep subject and region labels concise. The report hierarchy adds folders and
filenames to this root; see [Windows output paths](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#windows-output-paths).

The reports use the same structure as the bundled example, with results calculated from your inputs. Files will also be generated in your output directory: reports, connection density NIfTI files, figures, and tabular-formatted data. Note, this analysis assigns targets from the full streamline path and uses the fNode mask itself as a destination for `local` connection assignment. `local` and `local_search` mean different things. `local` records whether the entire streamline stays inside the fNode; `local_search` assigns a target from the two endpoints and a physical search radius. Changing that radius does not change the local-region mask.

Inspect the spatial previews and 3-D QC before interpreting results. FibrTraQ rejects incompatible image grids rather than silently resampling them. For `.trk`, when `axis_flip` is omitted, automatic mapping also stops if the best candidate flips any axis: sampled ROI overlap cannot establish anatomical alignment. For `.tck`, the mapping follows the format's world-coordinate convention without an orientation search. Every dataset still needs independent review. Follow the [spatial-verification steps](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/spatial-verification.md) to check the source geometry, pin the mapping and locate its metadata. Without `reference_image`, the grayscale background is tractogram-derived density, not independent anatomical evidence.

This analysis can be repeated for entire cohorts for subject-wise data.

### Optional settings

Add these arguments to `ftq_run_files()` when they apply:

| Argument | Use |
|---|---|
| `local_mask = file.path(node_dir, "local-region.nii.gz")` | Define a local region other than the focal mask |
| `endpoint_contexts = c(brain = file.path(node_dir, "brain-mask.nii.gz"))` | Check endpoint containment in trusted anatomy |
| `bvals = file.path(ref_dir, "diffusion.bval")` | Supply one b-value for every volume of a 4-D reference |
| `assignment_mode = "endpoint_voxels"` | Require an unambiguous focal/target endpoint-voxel pair |
| `assignment_mode = "local_search", search_radius_mm = 1` | Require an unambiguous endpoint pair within the specified physical radius; choose the radius for your study |
| `track_space = "voxelmm", axis_flip = c(FALSE, FALSE, FALSE)` | Pin a verified `.trk` mapping; this example requires confirmed voxel-millimeter coordinates with no axis flips. Native `.tck` uses `"world"` |
| `report_slices = c(20L, 30L, 40L)` | Select existing one-based axial slices for density heatmap previews; automatic selection is the default |

The [interpretation guide](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/interpretation.md#target-assignment-modes)
defines the assignment rules and their limits. `overwrite = FALSE` protects
existing outputs by default; choose a new output directory for another run.

<a id="what-fibrtraq-writes"></a>
<a id="interpret-results-correctly"></a>

## Read the main results

The main report presents the **RSC summary**, core maps, complete voxel data and
fate/coverage summaries first. Detailed endpoint QC, target patterns,
validation and provenance remain available in supporting sections.

**RSC measures focal-voxel coverage for one target:**

```text
RSC (%) = 100 × focal voxels with at least one association to this target
               / all focal voxels
```

The denominator includes focal voxels with zero assigned streamlines. One
focal voxel can be associated with several targets, so **target RSC values do
not have to sum to 100%**. The `local` connection row can overlap target rows.
RSC is not the percentage of streamlines assigned to a target.

For example, these are the two target rows from the bundled synthetic fixture,
which has two focal voxels:

| Target | Unique streamline connections | Connected focal voxels | RSC (%) |
|---|---:|---:|---:|
| cortex | 6 | 2 | 100 |
| brainstem | 4 | 2 | 100 |

Both RSC values are 100% because both focal voxels have associations with both
targets. The unique-streamline counts answer a different question and differ
between targets. These values come from the bundled [known-answer tables](inst/extdata/fibrtraq-known-answer/expected).

| What you want | Where to look |
|---|---|
| Unique streamline counts and target coverage | Main RSC summary; **Destination** identifies each target or the local comparison |
| Where associations occur inside the focal ROI | Core NIfTI maps and spatial previews |
| Every focal voxel and its counts | Complete per-voxel CSV; the readable all-voxel export adds coordinates and a column-role companion |
| All focal voxels ordered by a target or metric | Ranked CSV downloads, including zero-count voxels |
| Primary target, diversity and entropy patterns | Target report |
| Local, nonlocal target-positive and nonlocal target-negative streamlines | Unique-streamline fate summary; read its stated denominator |
| Parameter values, input provenance and computation checks | Readable Methods and Validation sections, with JSON downloads |
| Exact field definitions, units, aliases and missing values | The report's output guide, also [included with the package](inst/output-guide/OUTPUTS_README.md) |

A streamline can contribute to multiple focal voxels. **Do not sum per-voxel
counts to obtain unique-streamline totals.** Local fraction is `local / tracks`
where `tracks > 0`; zero-assignment voxels are unevaluable and have a separate
map validity mask. Ranked CSVs contain all focal voxels, rather than a top-50
subset. Use the original `*_per-voxel-hits.csv` for general
analysis; the coordinate-enriched view is for inspection and coordinate lookup.

Open `fibrtraq-report-index.html` first. Within each subject/focal-region
folder, the Target report and its five images are in `figures/`; its
`*_fibrtraq_exploratory-primary-target-lut.csv` lookup table is alongside that
folder. The main report remains under `reports/`, with images in
`report-figures/` and complete ranked CSVs in `tables/`. Core maps and any
tractogram-density background are in `maps/`. A subject's `prepared/` folder
exists only when the workflow needs to save a derived image or mask.

The standard workflow omits the extra exploratory CSVs and NIfTI maps.
Their calculations and explicit research exports remain available through
`ftq_exploratory_summary()`; see the [output inventory](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#what-fibrtraq-writes).

Computational validation, anatomical QC and report-generation status are
separate. A computational pass is not proof of biological validity. The 3-D
viewer shows a quota-limited subset for inspection; use full-run tables and
maps for quantitative claims. See [the full interpretation guide](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/interpretation.md)
and [output inventory](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#what-fibrtraq-writes).

## Run a regional study

Use `ftq_run_regional_study()` for multiple subjects with the same named set
of regions. Each requested region becomes focal in turn and the remaining
regions become targets. A two-subject, four-region study therefore creates
eight runs but still represents two biological subjects.

The [complete regional-study example](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#complete-regional-study-example)
shows the named input lists and the exact function call. Its default
`local_masks = "union"` defines local using the union of the supplied regions;
this differs from the one-analysis template's focal-mask default.

<a id="manifest-workflow"></a>
<a id="core-counting-function"></a>

For other workflows:

- [Manifest workflow](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#manifest-workflow): save a YAML/JSON configuration that can be shared and rerun.
- [Core in-memory function](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#core-counting-function): work directly with arrays and tractography objects.
- [S1c mouse-study recipe](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/mouse-study-workflow.md): use the verified filenames and audited settings for that dataset.

<a id="sharing-and-opening-reports"></a>

## Share the reports

Copy or ZIP the **whole results folder**, including subfolders. Extract it
before opening `fibrtraq-report-index.html`. Relative HTML links allow the
complete folder to move between Windows, macOS and Linux. Original source
paths remain in CSV/JSON provenance, so review identifiers and paths before
sharing a bundle publicly.

The report's output guide explains every table family and map. Older or
customized guides are retained when a newer package creates a current guide.
For historical Windows summaries containing `%5C` in broken links, use the
[summary-repair instructions](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#sharing-and-opening-reports).

## Common problems

| Problem | What to check |
|---|---|
| GitHub installation returns HTTP 404 | Confirm repository access and an authorized GitHub credential |
| Functions or arguments are missing after installation | Restart R; check the [loaded package revision](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#verify-an-installed-version) and reinstall the intended version |
| No focal hits | Review the tractogram coordinate convention, focal mask and selected mapping |
| Automatic `.trk` mapping stops because its best candidate flips an axis | Verify the tractogram-to-mask alignment independently, then explicitly set `track_space` and `axis_flip`; a higher ROI-overlap score does not prove the orientation is correct |
| A `.tck` run rejects geometry or coordinate settings | Use `track_space = "auto"` or `"world"` without extra flips, and resolve missing units or invalid/conflicting NIfTI transforms at their source; see [TCK requirements](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#mrtrix-tck-inputs) |
| Images have the same dimensions but are rejected | Check their NIfTI affines with `ftq_mask_grid_report()` |
| A 4-D image is rejected | Supply a b-value for every image volume |
| Existing-output error | Choose a new directory or deliberately set `overwrite = TRUE` |
| A Windows run stops because an output path is too long | Choose a shorter output directory and concise subject/region labels. The preflight error identifies the path that exceeds the supported length; see [Windows output paths](https://github.com/castrattonDSU/FibrTraQ/blob/main/docs/advanced-workflows.md#windows-output-paths) |
| A Windows report opens but a figure or download does not | Move the complete output folder to a short location such as `C:/FibrTraQ/`. Some Windows applications retain a [260-character full-path limit](https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation), even when R created the file successfully. Keep the report folder and its contents together |
| A report or viewer fails after the analysis completes | Open the saved run index: completed scientific files remain available and the failed stage is recorded |
| Anatomical QC is unavailable | Supply a trusted same-grid endpoint context mask; image-bound coverage alone is not anatomical evidence |
| Streamlines appear outside anatomy | Review complete paths, coordinate mapping and the input tractography; the viewer does not clip paths to a mask |

Open `?ftq_run_files`, `?ftq_run_regional_study` or `help(package = "FibrTraQ")`
for argument details. For support, open a [GitHub issue](https://github.com/castrattonDSU/FibrTraQ/issues)
with your package version, operating system, exact error and a minimal example.

## Validation, citation and license

The bundled fixture verifies reproducible software behavior. Unit and
integration tests cover counting semantics, coordinate mapping, streamed
readers, endpoint classes and R/C++ parity. GitHub Actions runs package checks
on Windows, macOS and Linux. Completed evidence and remaining acceptance gates
are recorded in [VALIDATION.md](https://github.com/castrattonDSU/FibrTraQ/blob/main/VALIDATION.md); larger studies are distinct from
the bundled synthetic example.

Use `citation("FibrTraQ")` or [CITATION.cff](https://github.com/castrattonDSU/FibrTraQ/blob/main/CITATION.cff) for citation metadata.
The associated manuscript is in preparation. FibrTraQ and its package-authored
synthetic fixture are licensed under the
[GNU General Public License, version 2 or later](https://github.com/castrattonDSU/FibrTraQ/blob/main/LICENSE.md).
Third-party dependencies retain their own licenses. Earlier copies distributed
under MIT, including release candidate `v0.3.0-rc2`, retain those terms.
See [CONTRIBUTING.md](https://github.com/castrattonDSU/FibrTraQ/blob/main/CONTRIBUTING.md) and [SUPPORT.md](https://github.com/castrattonDSU/FibrTraQ/blob/main/SUPPORT.md) for
contributor and support guidance.
