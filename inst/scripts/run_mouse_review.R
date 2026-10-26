# S1c mouse report review: bundled example -> one analysis -> Results.
# Install from GitHub and restart R first; see docs/mouse-study-workflow.md.
# Open a COPY of this script in RStudio and edit only these two paths.
# On Windows use forward slashes, for example "D:/FibrTraQ-data".
data_root <- "/path/to/FibrTraQ-data"
output_dir <- "/path/to/FibrTraQ-results/S1c-report-review-01"

library(FibrTraQ)

# 1. Locate the existing S1c files and protect earlier results ----
data_root <- path.expand(data_root)
output_dir <- path.expand(output_dir)
subject_id <- "S1c"

mask_dir <- file.path(data_root, "harmonized-data", "nodes")
focal_mask <- file.path(mask_dir, "S1c_cerebellum_native_ortho.nii.gz")
targets <- c(
  cerebrum = file.path(mask_dir, "S1c_cerebrum_native_ortho.nii.gz"),
  pontine = file.path(mask_dir, "S1c_pontine_native_ortho.nii.gz"),
  medullaspinal = file.path(mask_dir, "S1c_medullaspinal_native_ortho.nii.gz")
)
tractogram <- file.path(data_root, "real_validation", "trk", "S1c", "S1__B.trk")
reference_image <- file.path(
  data_root, "real_validation", "T2_nifti_raw",
  "S1c_m122_KH_JianLi_121120_63_1_1_f32.nii.gz"
)

required_files <- c(tractogram, focal_mask, targets, reference_image)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Check data_root; these input files were not found:\n",
       paste(missing_files, collapse = "\n"), call. = FALSE)
}
if (file.exists(output_dir)) {
  stop("Choose a NEW output_dir; this path already exists:\n",
       output_dir, call. = FALSE)
}

# 2. Check the installation using the README's bundled example ----
demo <- ftq_run_example(
  output_dir = file.path(tempdir(), "fibrtraq-known-answer"),
  overwrite = TRUE
)
stopifnot(demo$validation$passed)

# 3. Run S1c cerebellum against the other three regions ----
# The union defines the mouse study's local region: all four regions together.
local <- ftq_make_union_mask(
  masks = c(cerebellum = focal_mask, targets),
  template_image = focal_mask,
  output_path = file.path(output_dir, "S1c_regional-union.nii.gz")
)

analysis <- ftq_run_files(
  tractogram = tractogram,
  focal_mask = focal_mask,
  target_masks = targets,
  reference_image = reference_image,
  local_mask = local$path,
  output_dir = output_dir,
  subject = subject_id,
  focal_label = "cerebellum",
  assignment_mode = "local_search",
  search_radius_mm = 1,
  # Audited settings for these exact S1c files, not general defaults.
  track_space = "voxelmm",
  discretization = "round",
  axis_flip = c(FALSE, FALSE, FALSE),
  streamline_voxelization_mode = "segment_interpolation",
  # The independent dataset audit established voxel-index correspondence.
  # This copies the reference onto the mask header without resampling.
  reheader_image_to_mask_grid = TRUE
)

# 4. Check completion and open all reports from Results ----
stopifnot(
  all(analysis$index$run_status == "completed"),
  all(analysis$index$validation_passed),
  all(analysis$index$primary_report_status == "completed"),
  all(analysis$index$exploratory_report_status == "completed"),
  all(analysis$index$interactive_3d_qc_status == "completed")
)

print(analysis$index[, c(
  "primary_report", "exploratory_report", "interactive_3d_qc",
  "result_csv", "roi_summary_csv", "endpoint_summary_csv"
)])
cat("Mouse report index: ", analysis$report_index,
    "\nBundled example index: ", demo$report_index, "\n", sep = "")
browseURL(analysis$report_index)
# Choose Results. The main report links the Target report, 3-D viewer,
# complete tables, maps, and output definitions. Keep the whole output folder.
# The exploratory_report index field now points to the Target report.
# Ranked CSVs include all focal voxels, including zero counts.
# To inspect the synthetic example separately: browseURL(demo$report_index)
