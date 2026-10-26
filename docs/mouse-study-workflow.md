# Review one mouse brain: S1c cerebellum

This recipe runs the same public workflow as the [README](../README.md):
install FibrTraQ from GitHub, check the bundled example, run one analysis with
`ftq_run_files()`, and open the report index. It uses the **exact filenames and
audited settings for the study's S1c mouse**. These study inputs are not bundled
with the package; you need an authorized copy of the `FibrTraQ-data` folder.

## 1. Install from GitHub

In a fresh R session, run the README installation commands:

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

If the repository is private, your GitHub account needs access and R needs an
already configured GitHub credential. Use your existing authentication; do not
put an authentication key in the script. For a reproducible review, replace
`"main"` with the approved release tag or full commit SHA supplied with that
review. **Restart R after installation.**

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

Run the saved script from top to bottom. It first runs the small synthetic
known-answer example and stops if that validation fails. It then creates the
regional union mask and runs S1c with cerebellum as the focal ROI and cerebrum,
pontine, and medullaspinal as targets. All standard reports, maps, and tables
are generated.

The settings preserve the audited mouse workflow:

| Setting | Meaning for this study |
|---|---|
| Union of all four regional masks | Defines the `local` region; it is larger than the focal cerebellum |
| `assignment_mode = "local_search"`, `search_radius_mm = 1` | Uses the 1 mm endpoint search to assign target connections |
| `track_space = "voxelmm"`, `discretization = "round"`, no axis flips | Retains the established mapping for these tractogram and mask files |
| `streamline_voxelization_mode = "segment_interpolation"` | Includes voxels traversed between stored streamline points |
| `reheader_image_to_mask_grid = TRUE` | Copies this reference image onto the mask header, without resampling, after the dataset audit established voxel-index correspondence |

The reheader setting is specific to these audited files. For another dataset,
leave it `FALSE` unless an independent geometry audit supports that decision;
matching image dimensions alone is insufficient. The regional union is the
study's local classification mask, not a separately validated complete
brain-and-cord context mask. Computational validation does not establish
biological accuracy.

## 4. Review the reports

The script opens `fibrtraq-report-index.html` inside your chosen output folder.
Choose **Results** for `S1c / cerebellum`. From that page you can review:

- the RSC summary, spatial previews, complete tables, and maps;
- the **Target report** for supporting patterns and distributions;
- the interactive 3-D QC viewer for the sampled streamlines and ROI context;
- the output guide for column definitions and denominators.

The script checks that the analysis and all three reports completed, then
prints their paths. The synthetic example has its own index, also printed by
the script; open it with `browseURL(demo$report_index)` if desired.

The Target report and its five images are in `S1c/cerebellum/figures/`.
Its primary-target lookup CSV is in `S1c/cerebellum/`. The main report is
under `S1c/cerebellum/reports/`, with `report-figures/` for its images and
`tables/` for complete ranked CSVs. Every ranked CSV includes all focal
voxels, including those with zero counts. This mouse recipe needs a prepared
reference image for the audited reheader step, so a nonempty `S1c/prepared/`
folder is expected. Existing result folders are not changed.

Target RSC percentages describe overlapping focal-voxel memberships and do
not need to add to 100%. Read the [interpretation guide](interpretation.md)
alongside the report when comparing those values.

To share or move the results, keep the **entire output folder** together and
open its report index after copying or extracting it. Sending only the HTML
summary leaves its linked reports and figures behind. For another review run,
choose a new `output_dir`; the script preserves previous results.
