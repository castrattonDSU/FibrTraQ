# Verify spatial alignment before interpreting results

Use this procedure for each new dataset and after any change to registration,
tractography export, masks or coordinate settings. A completed run and a
computational-validation pass establish consistency of the calculated outputs;
they do not establish anatomical alignment.

The snippets below use the input variables from the
[single-subject README example](../README.md#run-fibrtraq-on-a-single-subject).
Use a new output directory for each review so earlier runs remain available.

## 1. Establish the source geometry

Confirm that the tractogram, focal mask, target masks and reference image belong
to the same subject and acquisition space. Retain the tractography-export
software/version, its coordinate convention and reference image, and the
registration transforms used to bring masks into that space. A matching filename
or array size does not establish this correspondence.

Use an anatomical reference independent of the streamline counts: for example,
the subject's diffusion b0 or a T1/T2 image registered to the diffusion grid.
The registration and image headers must themselves have documented provenance.
If the export convention or registration cannot be established, pause
quantitative interpretation and resolve that gap with the data producer.

## 2. Check the image grids

`ftq_mask_grid_report()` reads the images without changing them. This checks the
focal and target masks against the reference used in the README template:

```r
focal_path <- file.path(node_dir, "cerebellum.nii.gz")
reference_path <- file.path(ref_dir, "reference.nii.gz")
mask_paths <- c(focal = focal_path, targets)

grid_checks <- do.call(rbind, lapply(mask_paths, function(path) {
  ftq_mask_grid_report(mask = path, template_image = reference_path)
}))

grid_checks[, c("mask", "dims_match", "voxel_size_match",
                 "affine_match", "grid_match")]
stopifnot(all(grid_checks$grid_match), all(grid_checks$voxel_size_match))
```

Check any separate local, background or endpoint-context image too. For a
four-dimensional reference, first identify the appropriate three-dimensional
anatomical volume; see the [reference-image options](advanced-workflows.md#recommended-or-optional-files)
and `?ftq_run_files` for the package's b-value requirements.

**A grid match means that the read image geometry agrees.** It cannot detect
incorrectly copied headers, establish the subject's identity or prove that a
mask labels the intended anatomy. Stop on a mismatch or untrusted/missing spatial
metadata. Resolve registration at its source; do not copy headers or enable
`reheader_image_to_mask_grid` simply to make a check pass. That option is only
appropriate after independent proof of voxel-index correspondence.

## 3. Inspect independent anatomy and tractography

In a neuroimaging viewer that supports the exporter's coordinate convention,
overlay the focal and target masks on the trusted reference, then inspect the
tractogram in that same space. Record any transform the viewer applies.

- Inspect axial, coronal and sagittal planes at several levels, including the
  focal region and its boundaries. Check anterior/posterior, superior/inferior
  and left/right against acquisition geometry and identifiable landmarks.
- Use known asymmetric landmarks or acquisition-side information when
  available. Symmetric anatomy alone cannot settle a left/right reflection.
- Check for mirrored placement, shifts, scale differences, mask labels on the
  wrong structures, and unexplained streamline displacement. Full streamlines
  can extend outside a regional mask; containment in that mask is not the test.

Stop if the sources disagree or two coordinate conventions remain
indistinguishable at the resolution relevant to your analysis. More focal-mask
hits, plausible counts, or a visually convincing density image cannot choose
between them. Keep that uncertainty explicit rather than choosing the setting
that produces the preferred result.

## 4. Pin the supported mapping in the analysis

For `.trk`, after the source geometry and independent review support a mapping, set both
`track_space` and `axis_flip` explicitly in `ftq_run_files()`. Supported spaces
are `"voxel"`, `"voxelmm"`, `"vox10"`, `"world"` and `"trackvis"`; consult
`?FibrTraQ` for their definitions. Axis-flip entries correspond to the three
voxel-array axes and are not interchangeable with a viewer's display labels.

For example, **only when your evidence establishes voxel-millimeter coordinates
with no axis reflections**, use this replacement for the README's analysis
call. It keeps the README's assignment and discretization settings:

```r
analysis <- ftq_run_files(
  tractogram = file.path(trk_dir, "tractogram.trk"),
  focal_mask = focal_path,
  target_masks = targets,
  reference_image = reference_path,
  output_dir = output_dir,
  subject = subject_id,
  focal_label = "cerebellum",
  assignment_mode = "full_streamline",
  discretization = "any",
  track_space = "voxelmm",
  axis_flip = c(FALSE, FALSE, FALSE)
)
```

Record and pin `discretization` and the assignment/voxelization settings; they
affect counts independently of orientation. Retain the settings justified by your
analysis protocol rather than changing them to improve visual appearance.

For `.trk`, when flips are omitted, FibrTraQ stops before counting if its preferred
automatic candidate reflects an axis. The error identifies the candidate; it
does not certify it. Explicit reflections remain supported when independently
justified. An automatic unreflected choice or explicitly supplied wrong mapping
can still run, so **no flips is not a universal correction**.

For native `.tck`, use the same file workflow with the `.tck` filename and
`track_space = "auto"` or `"world"`. Leave `axis_flip = NULL` or set all three
entries to `FALSE`. TCK coordinates are RAS+ world millimeters; the mapping is
fixed by the NIfTI affine rather than selected by ROI overlap. Do not copy
voxel-millimeter or flip settings from a `.trk` analysis into a TCK run.

The reference NIfTI, or focal NIfTI when a reference is omitted, must declare
spatial units and a valid, invertible qform or sform. Declared qform and sform
transforms must agree when both are present. The run rejects missing units,
missing/invalid transforms and conflicting forms before counting. These
checks detect inconsistent metadata, but a consistent incorrect header still
requires the independent review above. Negative axes and oblique orientations
encoded by a valid affine are supported; no additional flips are needed.
See [TCK inputs](advanced-workflows.md#mrtrix-tck-inputs) for the supported
file encodings and geometry requirements.

## 5. Check the saved mapping and review the outputs

After the run, inspect the statuses and recorded selection:

```r
analysis$index[, c("subject", "focal_region", "analysis_status",
                   "validation_status", "selected_track_space",
                   "selected_axis_flip")]

metadata_path <- analysis$index$metadata_json[[1L]]
metadata_path
browseURL(analysis$report_index)
```

For `ftq_run_files()`, the index contains the single requested analysis. For a
regional study, review every row and use its corresponding `metadata_json`
path. The metadata JSON lives alongside the core CSVs in the subject/focal-region
folder and is linked from the report's Methods section. It records package
identity, input inventory, parameters, preflight checks and result summaries.
A run rejected before counting has a failed index entry and no scientific
run-metadata JSON; read its recorded error instead.

Confirm that the recorded mapping is the one you reviewed. Compare the report's
spatial previews and 3-D QC with the independent anatomical review. The 3-D
viewer shows a sampled subset, so use full-run maps and tables for quantitative
claims. Without a reference image, the grayscale report background is derived
from the same mapped streamlines and cannot independently confirm orientation.

Keep the input identifiers/checksums, software revision, mapping decision,
supporting review images, reviewer/date and any unresolved limitations with
the analysis. Repeat the review when inputs or transformations change, and
recompute affected outputs in a new folder if a mapping error is found.

See [validation evidence](../VALIDATION.md) for what the automated checks establish
and [coordinate safeguards](advanced-workflows.md#coordinate-and-grid-safeguards)
for additional grid-handling details.
