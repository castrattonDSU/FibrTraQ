# dTOR-985 left-thalamus FibrTraQ start-to-finish workflow.
#
# Run this from RStudio with:
# source(system.file("scripts", "fibrtraq-dtor985-start-to-finish.R", package = "FibrTraQ"))
#
# Set FIBRTRAQ_REPO_ROOT and FIBRTRAQ_DATA_ROOT if the repository and data are
# not adjacent directories.

# Set this to TRUE for the real workflow; set it to FALSE only when checking that paths and syntax are valid.
RUN_FULL_ANALYSIS <- TRUE

# Set this to TRUE when you want existing output files to be replaced by a fresh run.
OVERWRITE_OUTPUTS <- TRUE

# Set this to TRUE to open one results index after the workflow finishes.
OPEN_OUTPUTS <- TRUE

# Use "cpp" for the faster compiled backend. Change this to "r" only when you want the reference implementation.
COUNTING_BACKEND <- "cpp"

# Set this to TRUE so the long tractogram pass prints progress messages in the RStudio console.
SHOW_PROGRESS <- TRUE

# This controls how often the main tractogram pass prints progress when SHOW_PROGRESS is TRUE.
PROGRESS_EVERY <- 100000L

# This path is the local FibrTraQ source repository used for development installs and exact script locations.
repo_root <- Sys.getenv("FIBRTRAQ_REPO_ROOT", unset = normalizePath(getwd(), mustWork = FALSE))

# This path is the public-validation data root that contains the dTOR-985 raw data, ROIs, manifests, and outputs.
data_root <- Sys.getenv("FIBRTRAQ_DATA_ROOT", unset = file.path(dirname(repo_root), "FibrTraQ-data"))
dtor_root <- Sys.getenv("FIBRTRAQ_DTOR_ROOT", unset = file.path(data_root, "public-validation", "dtor985"))

# This path is the smaller pilot manifest; it tests the same workflow on a limited target set before the full run.
pilot_manifest <- file.path(dtor_root, "manifests", "dtor985-left-thalamus-harvardoxford-pilot.json")

# This path is the full manifest; it runs the left-thalamus ROI against all Harvard-Oxford target masks.
full_manifest <- file.path(dtor_root, "manifests", "dtor985-left-thalamus-harvardoxford-full.json")

# This is the raw TrackVis tractogram that FibrTraQ streams from disk during the pilot and full runs.
tractogram_path <- file.path(dtor_root, "raw", "dTOR_full_tractogram.trk")

# This is the focal ROI; every output row is one voxel from this left-thalamus mask.
focal_roi_path <- file.path(dtor_root, "rois", "harvard-oxford-thalamus", "focal", "Left_Thalamus.nii.gz")

# This is the anatomical background used for report overlays and 3D brain context.
background_image_path <- file.path(dtor_root, "anat", "MNI152_T1_2mm_on_HarvardOxford.nii.gz")

# This table lists the Harvard-Oxford target mask labels and their paths; it is useful for checking what was tested.
target_table_path <- file.path(dtor_root, "rois", "harvard-oxford-thalamus", "harvard-oxford-targets.csv")

# This directory stores the full FibrTraQ per-voxel CSV, RSC summary, fate summary, validation JSON, and metadata.
full_output_dir <- file.path(dtor_root, "outputs", "harvard-oxford-left-thalamus", "full")

# The target lookup sits beside the core results; the target report and images
# are written in this directory's figures/ subdirectory.
exploratory_output_dir <- full_output_dir

# This directory stores compact statistical outputs created for validation review.
statistics_output_dir <- file.path(dtor_root, "outputs", "harvard-oxford-left-thalamus", "full-statistics")

# This prefix ties target-report filenames to this dataset, ROI, atlas, and full run.
report_prefix <- "sub-dtor985_desc-left-thalamus-harvardoxford-full"

# These target labels are used only for the visual 3D QC sampling, not for limiting the quantitative full analysis.
qc_target_labels <- c(
  "Left_Frontal_Pole",
  "Left_Superior_Frontal_Gyrus",
  "Left_Precentral_Gyrus",
  "Left_Postcentral_Gyrus",
  "Left_Putamen",
  "Left_Brain_Stem"
)

# This helper prints section headers so the RStudio console is easier to follow during a long workflow.
announce <- function(text) {
  cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  cat(text, "\n")
  cat(paste(rep("=", 78), collapse = ""), "\n", sep = "")
}

# This helper stops early if an input path is missing; that prevents a long run from failing halfway through.
require_existing_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path, call. = FALSE)
  invisible(path)
}

# This helper creates output directories before writing CSVs, reports, and figures.
require_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# Use the package browser opener for portable local-file URLs.
open_local_file <- function(path) {
  FibrTraQ:::ftq_qc_open_output(path)
}
# This block loads the local development copy when possible, so RStudio uses the code in this repository.
announce("Step 1: Load FibrTraQ and required helper packages")
if (requireNamespace("pkgload", quietly = TRUE) && file.exists(file.path(repo_root, "DESCRIPTION"))) {
  pkgload::load_all(repo_root)
} else {
  library(FibrTraQ)
}
cat("Counting backend selected:", COUNTING_BACKEND, "\n")

# Installed companion scripts work without a source checkout. A source fallback
# also supports development copies that have not been installed yet.
streamline_figure_script <- system.file(
  "scripts", "fibrtraq-dtor985-plot-streamlines.R", package = "FibrTraQ"
)
if (!nzchar(streamline_figure_script)) {
  streamline_figure_script <- file.path(
    repo_root, "inst", "scripts", "fibrtraq-dtor985-plot-streamlines.R"
  )
}

# jsonlite is needed to read the validation JSON and summarize whether internal validation passed.
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite before running this script.", call. = FALSE)

# RNifti is needed by FibrTraQ for NIfTI reading/writing and by the report generation code.
if (!requireNamespace("RNifti", quietly = TRUE)) stop("Install RNifti before running this script.", call. = FALSE)

# ggplot2 is needed by the streamline figure script for static and anatomical overlay figures.
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 before running this script.", call. = FALSE)

# This block verifies the exact raw inputs before any analysis starts.
announce("Step 2: Verify raw inputs and script paths")
require_existing_file(pilot_manifest, "Pilot manifest")
require_existing_file(full_manifest, "Full manifest")
require_existing_file(tractogram_path, "Raw dTOR-985 tractogram")
require_existing_file(focal_roi_path, "Left-thalamus focal ROI")
require_existing_file(background_image_path, "MNI anatomical background image")
require_existing_file(target_table_path, "Harvard-Oxford target table")
require_existing_file(streamline_figure_script, "Streamline figure script")
require_directory(full_output_dir)
require_directory(exploratory_output_dir)
require_directory(statistics_output_dir)

# This reads the target table so the console records how many atlas targets are available.
target_table <- utils::read.csv(target_table_path, stringsAsFactors = FALSE)
cat("Target masks available:", nrow(target_table), "\n")

# This preflight checks coordinate settings, ROI grids, target masks, output collisions, and path validity for the pilot.
announce("Step 3: Preflight the pilot manifest")
pilot_preflight <- ftq_preflight(pilot_manifest, overwrite = OVERWRITE_OUTPUTS)
print(pilot_preflight$checks)
if (!isTRUE(pilot_preflight$passed)) stop("Pilot preflight failed. Fix the failed error checks above before running.", call. = FALSE)

# This preflight performs the same checks for the full 117-target run.
announce("Step 4: Preflight the full manifest")
full_preflight <- ftq_preflight(full_manifest, overwrite = OVERWRITE_OUTPUTS)
print(full_preflight$checks)
if (!isTRUE(full_preflight$passed)) stop("Full preflight failed. Fix the failed error checks above before running.", call. = FALSE)

# This switch lets the same file be used for path/syntax checks without launching the full tractogram workflow.
if (!isTRUE(RUN_FULL_ANALYSIS)) {
  announce("RUN_FULL_ANALYSIS is FALSE, so the script stopped after path and preflight checks")
  stop("Set RUN_FULL_ANALYSIS <- TRUE at the top of the script to run the full workflow.", call. = FALSE)
}

# The pilot run is a fast sanity check with a small target set before the full 117-target analysis.
announce("Step 5: Run the pilot FibrTraQ analysis")
pilot_run <- ftq_run(
  pilot_manifest,
  overwrite = OVERWRITE_OUTPUTS,
  backend = COUNTING_BACKEND,
  progress = SHOW_PROGRESS,
  progress_every = PROGRESS_EVERY,
  profile = TRUE
)
cat("Pilot per-voxel result CSV:\n", pilot_run$outputs$result_csv, "\n", sep = "")
cat("Pilot RSC summary CSV:\n", pilot_run$outputs$rsc_csv, "\n", sep = "")
cat("Pilot validation JSON:\n", pilot_run$outputs$validation_json, "\n", sep = "")

# The full run streams the raw .trk file and counts each focal voxel's connections to all atlas targets.
announce("Step 6: Run the full FibrTraQ analysis")
full_started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
full_run <- ftq_run(
  full_manifest,
  overwrite = OVERWRITE_OUTPUTS,
  backend = COUNTING_BACKEND,
  progress = SHOW_PROGRESS,
  progress_every = PROGRESS_EVERY,
  profile = TRUE
)
cat("Full per-voxel result CSV:\n", full_run$outputs$result_csv, "\n", sep = "")
cat("Full RSC summary CSV:\n", full_run$outputs$rsc_csv, "\n", sep = "")
cat("Full fate summary CSV:\n", full_run$outputs$fate_csv, "\n", sep = "")
cat("Full validation JSON:\n", full_run$outputs$validation_json, "\n", sep = "")

# This reads the validation JSON so the script fails immediately if FibrTraQ detected an internal validation problem.
full_validation <- jsonlite::fromJSON(full_run$outputs$validation_json, simplifyVector = FALSE)
if (!isTRUE(full_validation$passed)) stop("Full run validation JSON did not pass.", call. = FALSE)
cat("Full run validation passed:", full_validation$passed, "\n")

# This report summarizes the per-voxel table with target totals and spatial figures.
announce("Step 7: Build the target report and spatial figures")
exploratory_report <- ftq_exploratory_report(
  result_csv = full_run$outputs$result_csv,
  output_dir = exploratory_output_dir,
  template_image = full_run$manifest$fnode,
  background_image = full_run$manifest$background_image,
  prefix = report_prefix,
  assignment_mode = full_run$manifest$target_intersection_mode,
  search_radius_mm = full_run$manifest$local_search_radius_mm,
  overwrite = OVERWRITE_OUTPUTS
)
cat("Target HTML report:\n", exploratory_report$html, "\n", sep = "")
cat("Target-report figure paths:\n")
print(exploratory_report$figures)

# This standalone viewer samples real streamlines and renders rotatable 3D QC with the whole focal ROI volume.
announce("Step 8: Build the interactive 3D QC viewer from the full run")
qc_viewer <- ftq_interactive_qc_viewer(
  full_run,
  output_html = file.path(exploratory_output_dir, "figures", paste0(report_prefix, "_fig_interactive-3d-qc.html")),
  target_labels = qc_target_labels,
  sample_per_target = 12L,
  max_streamlines = 3000000L,
  max_points_per_streamline = 140L,
  brain_points = 12000L,
  overwrite = OVERWRITE_OUTPUTS,
  open = FALSE
)
cat("Interactive 3D QC viewer:\n", qc_viewer$html, "\n", sep = "")
print(qc_viewer$summary)

# This script creates the additional streamline evidence figures used in the validation report.
announce("Step 9: Generate representative streamline figures and anatomical overlay")

# This environment passes OPEN_OUTPUTS into the figure script without polluting unrelated RStudio variables.
figure_script_env <- new.env(parent = globalenv())

# The parent workflow opens one results index after all companion steps finish.
figure_script_env$OPEN_OUTPUTS <- FALSE

# This runs the validated figure script inside the controlled environment above.
source(streamline_figure_script, local = figure_script_env)

# This summary collapses the per-voxel table into ROI-level, target-level, and voxel-level statistics.
announce("Step 10: Write publication-oriented summary statistics")
subject_summary <- ftq_stat_subject_summary(
  result_csv = full_run$outputs$result_csv,
  subject = "sub-dtor985",
  group = "public-validation",
  roi = "Harvard-Oxford left thalamus"
)

# This estimates target-fraction gradients across focal-voxel x/y/z coordinates for exploratory spatial structure.
spatial_gradients <- ftq_stat_spatial_gradients(
  result_csv = full_run$outputs$result_csv,
  min_connected_voxels = 10L
)

# These output paths keep the statistical summaries separate from the main FibrTraQ result tables.
roi_summary_csv <- file.path(statistics_output_dir, paste0(report_prefix, "_roi-summary.csv"))
target_summary_csv <- file.path(statistics_output_dir, paste0(report_prefix, "_target-summary.csv"))
voxel_summary_csv <- file.path(statistics_output_dir, paste0(report_prefix, "_voxel-summary.csv"))
spatial_gradients_csv <- file.path(statistics_output_dir, paste0(report_prefix, "_spatial-gradients.csv"))

# These CSVs support independent review of the quantitative validation.
utils::write.csv(subject_summary$roi, roi_summary_csv, row.names = FALSE)
utils::write.csv(subject_summary$targets, target_summary_csv, row.names = FALSE)
utils::write.csv(subject_summary$voxels, voxel_summary_csv, row.names = FALSE)
utils::write.csv(spatial_gradients, spatial_gradients_csv, row.names = FALSE)

# These console prints show the main quantitative interpretation without opening any files.
cat("ROI summary CSV:\n", roi_summary_csv, "\n", sep = "")
cat("Target summary CSV:\n", target_summary_csv, "\n", sep = "")
cat("Voxel summary CSV:\n", voxel_summary_csv, "\n", sep = "")
cat("Spatial gradients CSV:\n", spatial_gradients_csv, "\n", sep = "")
cat("\nROI-level summary:\n")
print(subject_summary$roi)
cat("\nTop 20 target-level summaries:\n")
print(utils::head(subject_summary$targets, 20))
cat("\nTop 20 spatial gradients by smallest q-value:\n")
print(utils::head(spatial_gradients[order(spatial_gradients$q_value), ], 20))

# These are the key files expected at the end of a successful start-to-finish run.
announce("Step 11: Verify final outputs exist")
expected_outputs <- c(
  full_run$outputs$result_csv,
  full_run$outputs$rsc_csv,
  full_run$outputs$fate_csv,
  full_run$outputs$validation_json,
  exploratory_report$html,
  exploratory_report$summary$outputs$lut_csv,
  unlist(exploratory_report$figures, use.names = FALSE),
  qc_viewer$html,
  roi_summary_csv,
  target_summary_csv,
  voxel_summary_csv,
  spatial_gradients_csv,
  file.path(exploratory_output_dir, "figures", "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines.png"),
  file.path(exploratory_output_dir, "figures", "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines-anatomical-overlay.png"),
  file.path(exploratory_output_dir, "figures", "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines-3d-interactive.html")
)
missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
if (length(missing_outputs)) stop("Missing expected output(s):\n", paste(missing_outputs, collapse = "\n"), call. = FALSE)
cat("All expected outputs exist.\n")

# A single landing page keeps the complete reports and viewer discoverable.
results_dir <- dirname(full_output_dir)
results_index <- file.path(results_dir, "fibrtraq-report-index.html")
results_csv <- file.path(results_dir, "fibrtraq-run-index.csv")
FibrTraQ:::ftq_check_overwrite(c(results_index, results_csv), OVERWRITE_OUTPUTS)
results_row <- FibrTraQ:::ftq_regional_study_run_row(
  subject = "sub-dtor985", focal_region = "left thalamus", run = full_run,
  exploratory = exploratory_report, viewer = qc_viewer,
  started_at = full_started_at, error = NULL
)
utils::write.csv(results_row, results_csv, row.names = FALSE)
FibrTraQ:::ftq_write_regional_study_index(
  results_row, path = results_index, output_dir = results_dir,
  description = report_prefix
)
if (!is.null(full_run$report$html) && file.exists(full_run$report$html)) {
  FibrTraQ:::ftq_report_add_navigation(
    full_run$report$html,
    list("Results index" = results_index,
         "Target report" = exploratory_report$html,
         "3-D QC viewer" = qc_viewer$html)
  )
}
announce("Step 12: Open the results index")
if (isTRUE(OPEN_OUTPUTS)) {
  open_local_file(results_index)
} else {
  cat("OPEN_OUTPUTS is FALSE, so no browser windows were opened.\n")
}

# These final lines make the most important paths available after source() finishes in RStudio.
final_outputs <- list(
  report_index = results_index,
  per_voxel_csv = full_run$outputs$result_csv,
  rsc_summary_csv = full_run$outputs$rsc_csv,
  fate_summary_csv = full_run$outputs$fate_csv,
  validation_json = full_run$outputs$validation_json,
  exploratory_html = exploratory_report$html,
  interactive_qc_html = qc_viewer$html,
  roi_summary_csv = roi_summary_csv,
  target_summary_csv = target_summary_csv,
  voxel_summary_csv = voxel_summary_csv,
  spatial_gradients_csv = spatial_gradients_csv
)
cat("\nFinal output paths:\n")
print(final_outputs)

# Returning this list lets RStudio users inspect the paths with final_outputs$exploratory_html, etc.
final_outputs
