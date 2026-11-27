# Review one mouse brain: S1c cerebellum

This recipe runs the same public workflow as the [README](../README.md):
run one analysis with `fibrtraq_single()` and open its HTML report. It uses the **exact filenames and
audited settings for the study's S1c mouse**. These study inputs are not bundled
with the package; you need an authorized copy of the `FibrTraQ-data` folder.

## 1. Install the matching package revision

Follow the [installation steps](../README.md#install) for the two-command package revision, then restart R. Run the bundled example in the README if you want to verify your installation first.

## 2. Open the script and set two folder paths

Open a copy of [run_mouse_review.R](../inst/scripts/run_mouse_review.R) in
RStudio. The installed package also includes this script; locate it with:

```r
system.file("scripts", "run_mouse_review.R", package = "FibrTraQ")
```

Use **File → Open File** to open that path, then **Save As** to keep an editable
copy in your own working folder. Change only the two lines at the top:

| Setting | What goes there |
|---|---|
| `data_root` | The full path to your existing `FibrTraQ-data` folder, containing both `real_validation` and `harmonized-data` |
| `output_dir` | The full path to a **new, unused** folder for this review, such as a folder named `S1c-report-review-01` inside your results folder |

The paths in the distributed script are explicitly marked templates. Replace
them with your workstation's full paths. On Windows use forward slashes, for
example `"D:/FibrTraQ-data"` and
`"D:/FibrTraQ-results/S1c-report-review-01"`. Spaces in folder names are allowed.
Keep the filenames below unchanged:

```text
FibrTraQ-data/
├── harmonized-data/nodes/
│   ├── S1c_cerebellum_native_ortho.nii.gz
│   ├── S1c_cerebrum_native_ortho.nii.gz
│   ├── S1c_pontine_native_ortho.nii.gz
│   └── S1c_medullaspinal_native_ortho.nii.gz
└── real_validation/
    ├── trk/S1c/S1__B.trk
    └── T2_nifti_raw/S1c_m122_KH_JianLi_121120_63_1_1_f32.nii.gz
```

`S1__B.trk` has **two underscores**. The script labels the subject explicitly
as `S1c`; it does not infer the subject from that tractogram filename. It checks
these six inputs and refuses an existing output folder before writing results.

## 3. Run the script

Run the saved script from top to bottom. It runs `fibrtraq_single()` with the
regional union mask and the cerebellum as the focal ROI, with cerebrum,
pontine, and medullaspinal as targets. All standard reports, maps, and tables
are generated in one output folder.

The settings preserve the audited mouse workflow:

| Setting | Meaning for this study |
|---|---|
| `local_mask = "union"` | Defines the `local` region; it is larger than the focal cerebellum |
| `assignment_mode = "local_search"`, `search_radius_mm = 1` | Uses the 1 mm endpoint search to assign target connections |
| `track_space = "voxelmm"`, `discretization = "round"`, no axis flips | Retains the established mapping for these tractogram and mask files |
| `streamline_voxelization_mode = "segment_interpolation"` | Samples additional points along segments between stored streamline points |
| `reheader_image_to_mask_grid = TRUE` | Copies this reference image onto the mask header, without resampling, after the dataset audit established voxel-index correspondence |

The reheader setting is specific to these audited files. For another dataset,
leave it `FALSE` unless an independent geometry audit supports that decision;
matching image dimensions alone is insufficient. The regional union is the
study's local classification mask, not a separately validated complete
brain-and-cord context mask. Computational validation does not establish
biological accuracy.

## 4. Review the reports

In RStudio the script opens the report in `S1c/cerebellum/`. Review the RSC table and density previews first. Select a destination in the table to open its count map, or open the 3D viewer to inspect sampled streamlines. All images are in `figures/`, CSVs in `tables/`, and core
NIfTI images in `maps/`. This recipe prepares the regional union mask and a reference in
`S1c/prepared/` for its audited reheader step.

For a quick review, choose one destination and read its unique streamline count, connected focal-voxel count and RSC percentage. Open that destination’s NIfTI from the table. Targets can share focal voxels, so their RSC percentages may sum to more than 100%. [Interpretation guide](interpretation.md)

To share the run, zip the directory printed under **Share this folder** and have the recipient extract it before opening the HTML. For another review run,
choose a new `output_dir`; the script preserves previous results.

## Orientation for the six harmonized controls

For **S15, S16, S17, S1c, S22 and S23**, comparisons with the original diffusion
b0 and T2 images support voxel-millimeter coordinates with no axis reflection
for the audited tractograms and harmonized regional masks. Before applying
those settings to a copy, verify its tractogram and mask checksums against
the audited files. This matters especially for renamed tractograms: the
audited S15 file is `dti2.trk` and the audited S1c file is `S1__B.trk`.
Repeat the [spatial review](spatial-verification.md) for different exports or masks.

For a single control, add `track_space = "voxelmm"` and
`axis_flip = c(FALSE, FALSE, FALSE)` to the existing `fibrtraq_single()` call.
Use a fresh output folder. These two arguments set only the coordinate mapping.
Keep the full-streamline assignment and stored-point voxelization defaults for
the standard analysis. The S1c recipe above uses a separate endpoint-search protocol.

For all six together, follow the [README study sheet example](../README.md#run-a-study)
and put only those six verified subjects in the sheet. With `sheet` and a fresh
`results_dir` set as in that example, use the same mapping for the check and run:

```r
checked <- fibrtraq_study(
  study_sheet = sheet, focal_regions = "cerebellum", output_dir = results_dir,
  track_space = "voxelmm", axis_flip = c(FALSE, FALSE, FALSE), check_only = TRUE
)
print(checked$index)
stopifnot(checked$passed)

study <- fibrtraq_study(
  study_sheet = sheet, focal_regions = "cerebellum", output_dir = results_dir,
  track_space = "voxelmm", axis_flip = c(FALSE, FALSE, FALSE), open = TRUE
)
```
