# FibrTraQ Extended Validation And Publication Evidence

This directory contains the developer-facing machinery for testing FibrTraQ
beyond its bundled known-answer example. It is not the package's getting-started
guide, and the presence of a dataset template does not mean that the associated
analysis has been completed.

The repository-level
[`VALIDATION.md`](https://github.com/castrattonDSU/FibrTraQ/blob/main/VALIDATION.md)
summarizes current release evidence. This directory provides the longer
protocols, manifest templates, scripts, and registry needed to extend that
evidence. It separates:

1. package-owned known-answer checks that require no download;
2. local regression checks against accepted, versioned outputs;
3. production-scale implementation benchmarks;
4. independent public benchmark applications that require external data and
   study-specific preparation.

See `publication-statistical-testing-plan.md` for the validation statistics,
negative controls, and biological hypothesis tests recommended for publication.
See `endpoint-analysis-data-dictionary.md` for the endpoint-versus-intersection,
target-penetration, reference-grid, and optional anatomical-context definitions.

## Evidence Bundle Policy

The R package should install only the small package-authored known-answer
fixture and lightweight validation assets: manifest templates, setup/run
scripts, dataset registry notes, benchmark summaries, and statistical testing
plans. External raw tractograms, downloaded atlases, derived masks, full HTML
reports, rendered figures, CSV outputs, and manuscript handoff files belong in
an external evidence bundle such as OSF, Zenodo, or a companion repository.

Each evidence bundle should include:

- a plain-language `README` describing the datasets, software version, and run
  order,
- copied manifests used for the final analyses,
- CSV summaries and validation JSON files,
- full HTML reports and PNG figures,
- checksums for raw and derived files when practical,
- notes identifying which analyses are ground-truth validation, negative
  controls, biological plausibility checks, and exploratory demonstrations.

The installed scripts should use configurable data roots or copied manifests
rather than private absolute paths whenever possible.

## Evidence Status

| Status | Dataset or check | Current evidence | Intended role |
| --- | --- | --- | --- |
| Complete | Bundled known-answer example | Seven exact scientific tables, report generation, integrity inventory, and automated tests | Primary user-runnable software-reproducibility check |
| Complete, historical scope | dTOR-985 normative connectome | Production-scale R/C++ comparison for per-voxel, RSC, and fate outputs on 11,820,000 streamlines | Backend equivalence and performance evidence for established core outputs |
| Template only | ISMRM 2015 Tractography Challenge | Manifest-generation and suite-runner support; no final release evidence is claimed here | Controlled tractography benchmark |
| Registry only | TractoInferno / OpenNeuro ds003900 | Candidate source and setup notes only | Real human, multi-site generalization |
| Template only | FiberCup phantom | Candidate manifest and setup notes only | Independent controlled phantom |
| Template only | Allen CCF plus compatible mouse dMRI tractography | Candidate manifest and setup notes only | Atlas-driven mouse application with external biological context |
| Registry only | HCP diffusion | Candidate source and setup notes only | High-quality human stress test |

`public-datasets.yml` is the machine-readable status registry. An item moves to
completed evidence only after the final run, artifacts, provenance, acceptance
checks, and limitations have been recorded. Templates and candidates must not
be cited as completed validation.

## Acceptance Criteria

Each manifest in a publication validation run should meet these criteria:

- `ftq_preflight()` passes with zero blocking errors.
- `ftq_run()` completes from raw tractogram and ROI files, not from summarized
  CSVs.
- `ftq_validate_run()` passes. If known-good outputs are listed in `expected`,
  per-voxel tables, RSC summaries, and maps match within the configured
  tolerance.
- The output table includes `leave_no_target` and all counts are non-negative.
- When endpoint analysis is enabled, endpoint classes partition every existing
  full-path target intersection, penetration histograms reconcile exactly to
  those intersections, and percentages identify their denominators in their
  column names.
- Endpoint analysis reports exact membership as the primary estimate and the
  fixed 26-neighbor, one-voxel dilation only as voxel-grid sensitivity.
- Computational validation passes independently of the biological-QC status.
  A biological status of `review` or `unavailable` must be resolved or disclosed
  before making anatomical-accuracy claims.
- The HTML report is generated and includes spatial figures. Anatomical
  background overlays are recommended for final figures, but they are not
  required for the software to run.
- Re-running the same manifest with the same inputs gives identical per-voxel,
  RSC, endpoint-summary, endpoint-penetration, and endpoint-context outputs.

## Running The Suite

From RStudio:

```r
library(FibrTraQ)

manifest_dir <- system.file("validation", "manifests", package = "FibrTraQ")
templates <- list.files(manifest_dir, full.names = TRUE)
templates

# After copying/editing manifests so the paths point to downloaded data:
summary <- ftq_run_validation_suite(
  manifests = "/path/to/edited-validation-manifests",
  output_csv = "/path/to/fibrtraq-validation-summary.csv",
  overwrite = TRUE
)
summary
```

From the shell:

```bash
Rscript "$(Rscript -e 'cat(system.file("scripts", "fibrtraq-validation-suite.R", package = "FibrTraQ"))')" \
  /path/to/edited-validation-manifests \
  --output-csv /path/to/fibrtraq-validation-summary.csv \
  --overwrite
```

Use `--preflight-only` first when setting up a new dataset.

## ISMRM 2015 Ground-Truth Manifests

When the ISMRM scoring data are unpacked, use the segmentation config to build
one manifest per configured bundle/sub-bundle:

```r
library(FibrTraQ)

data_root <- Sys.getenv("FIBRTRAQ_DATA_ROOT", unset = "/path/to/FibrTraQ-data")
scoring <- file.path(data_root, "public-validation", "ismrm2015", "scoring")
manifest_dir <- file.path(data_root, "public-validation", "ismrm2015", "manifests", "ground-truth")
output_root <- file.path(data_root, "public-validation", "ismrm2015", "outputs", "ground-truth")

generated <- ftq_make_ismrm2015_manifests(
  scoring_dir = scoring,
  output_dir = manifest_dir,
  run_output_root = output_root,
  force = TRUE
)

preflight <- ftq_run_validation_suite(manifest_dir, preflight_only = TRUE)
results <- ftq_run_validation_suite(
  manifest_dir,
  output_csv = file.path(output_root, "ismrm2015-ground-truth-validation-summary.csv"),
  overwrite = TRUE
)
```

Each generated manifest uses the public ground-truth `.trk` file for one bundle,
the configured head endpoint as the focal ROI, the configured tail endpoint as
the target ROI, and `t1.nii.gz` as the anatomical report background when present.

## dTOR-985 Exploratory Utility Test

dTOR-985 is the recommended high-impact utility demonstration because it is a
large normative, whole-brain MNI-space tractogram. Instead of asking whether one
known bundle reaches one known endpoint, this test asks a more realistic
exploratory question: from a focal left thalamus ROI, which cortical and
subcortical Harvard-Oxford targets are most represented at each focal voxel,
where connectivity is diverse, and where streamlines leave the focal ROI
without reaching an atlas target?

Prepare the Figshare metadata/download target:

```bash
Rscript "$(Rscript -e 'cat(system.file("scripts", "fibrtraq-dtor985-download.R", package = "FibrTraQ"))')" \
  /path/to/FibrTraQ-data/public-validation/dtor985 \
  --download
```

Prepare the Harvard-Oxford masks and manifests with a Python environment that
has `nilearn` and `nibabel` installed:

```bash
python "$(Rscript -e 'cat(system.file("scripts", "fibrtraq-dtor985-prepare.py", package = "FibrTraQ"))')" \
  --data-root /path/to/FibrTraQ-data \
  --force
```

Then run preflight and a pilot from RStudio:

```r
library(FibrTraQ)

data_root <- Sys.getenv("FIBRTRAQ_DATA_ROOT", unset = "/path/to/FibrTraQ-data")
dtor_root <- file.path(data_root, "public-validation", "dtor985")
pilot_manifest <- file.path(dtor_root, "manifests", "dtor985-left-thalamus-harvardoxford-pilot.json")
full_manifest <- file.path(dtor_root, "manifests", "dtor985-left-thalamus-harvardoxford-full.json")

ftq_preflight(pilot_manifest)
pilot <- ftq_run(pilot_manifest, overwrite = TRUE)

ftq_preflight(full_manifest)
full <- ftq_run(full_manifest, overwrite = TRUE)

explore <- ftq_exploratory_report(
  result_csv = full$outputs$result_csv,
  output_dir = file.path(dtor_root, "outputs", "harvard-oxford-left-thalamus", "exploratory-report"),
  template_image = full$manifest$fnode,
  background_image = full$manifest$background_image,
  prefix = "sub-dtor985_desc-left-thalamus-harvardoxford"
)
explore$html
```

The pilot manifest is for coordinate/report QC. The full manifest is the exact
many-target analysis; the exploratory report converts that full result into
primary-target, target-diversity, entropy, total-target-hit, and
leave-no-target maps.

The dTOR `.trk` stores raw TrackVis points, not ROI-grid voxel indices. The
generated manifests therefore use `track_space = "trackvis"` so FibrTraQ applies
the `.trk` header affine before converting streamlines to Harvard-Oxford voxels.
They also use `execution.track_reader = "stream"` so the full tractogram is read
in chunks rather than materialized in memory. No additional `axis_flip` should be
applied after this TrackVis transform.

## Dataset Registry

`public-datasets.yml` records each dataset's evidence status, source URLs,
intended role, and practical setup notes. The manifests in `manifests/` are
templates; they are not expected to pass until copied and edited to point at
local downloaded data, and they are not completed evidence by themselves.
