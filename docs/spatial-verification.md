# Verify spatial alignment before interpreting results

Use this procedure for each new dataset and after any change to registration,
tractography export, masks or coordinate settings. A completed run and a
computational-validation pass establish consistency of the calculated outputs;
they do not establish anatomical alignment.

The single-subject example below uses `subject_dir` from the
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

Check that the focal, target and local masks, reference image, background and endpoint-context masks are in the same subject space. Inspect their voxel dimensions, sizes and NIfTI transforms in your imaging software. Both public workflows check input geometry before counting and record their findings in the report’s **Preflight** section.

For a four-dimensional diffusion reference, supply one b-value for each volume so FibrTraQ can prepare a three-dimensional low-b image. See the [reference-image options](advanced-workflows.md#recommended-or-optional-files).

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
`track_space` and `axis_flip` explicitly in `fibrtraq_single()` or `fibrtraq_study()`. Supported spaces
are `"voxel"`, `"voxelmm"`, `"vox10"`, `"world"` and `"trackvis"`; consult
`?fibrtraq_single` for the mapping options. Axis-flip entries correspond to the three
voxel-array axes and are not interchangeable with a viewer's display labels.

For example, **only when your evidence establishes voxel-millimeter coordinates
with no axis reflections**, the single-subject README call becomes:

```r
analysis <- fibrtraq_single(
  tractogram = file.path(subject_dir, "tractogram.trk"),
  focal_mask = file.path(subject_dir, "cerebellum.nii.gz"),
  target_masks = c(cerebrum = file.path(subject_dir, "cerebrum.nii.gz"),
                   pontine = file.path(subject_dir, "pontine.nii.gz"),
                   medullaspinal = file.path(subject_dir, "medullaspinal.nii.gz")),
  reference_image = file.path(subject_dir, "reference.nii.gz"),
  output_dir = file.path(dirname(subject_dir), "001-verified-results"),
  subject = "001",
  focal_label = "cerebellum",
  track_space = "voxelmm",
  axis_flip = c(FALSE, FALSE, FALSE),
  open = TRUE
)
```

Record and pin `discretization` and the assignment/voxelization settings; they
affect counts independently of orientation. Retain the settings justified by your
analysis protocol rather than changing them to improve visual appearance.

With automatic TRK mapping, FibrTraQ scores candidate mappings on the first
200 streamlines by default, or all records in a smaller file. Their positions
and order can favor a reflection that puts more of the sample in the focal mask.
If the preferred candidate would reflect an axis and `axis_flip` was omitted,
the analysis stops before counting. The masks remain in their original orientation.
Version `0.3.5.9003` also reports this condition in `check_only = TRUE` for
each affected subject and focal region. Use independent anatomy to evaluate
the proposed mapping, even when the check passes. [Harmonized mouse controls](mouse-study-workflow.md#orientation-for-the-six-harmonized-controls)

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

For `fibrtraq_single()`, the index contains the single requested analysis. For a
regional study, review every row and use its corresponding `metadata_json`
path. The metadata JSON lives in the subject/focal-region
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
