# View the completed dTOR-985 left-thalamus FibrTraQ validation outputs.
#
# Run this from RStudio after installing/loading the package. It uses exact
# paths derived from FIBRTRAQ_DTOR_ROOT and opens one available results entry.
# Historical supporting files are optional when viewing an existing report.

data_root <- Sys.getenv("FIBRTRAQ_DATA_ROOT", unset = file.path(getwd(), "FibrTraQ-data"))
dtor_root <- Sys.getenv("FIBRTRAQ_DTOR_ROOT", unset = file.path(data_root, "public-validation", "dtor985"))
left_root <- file.path(dtor_root, "outputs", "harvard-oxford-left-thalamus")
full_dir <- file.path(left_root, "full")
legacy_report_dir <- file.path(left_root, "full-exploratory-report")
report_prefix <- "sub-dtor985_desc-left-thalamus-harvardoxford-full"
target_report <- file.path(full_dir, "figures", paste0(report_prefix, "_fibrtraq_target-report.html"))
legacy_report <- file.path(legacy_report_dir, paste0(report_prefix, "_fibrtraq_exploratory-report.html"))
# New runs use full/figures; historical reports remain readable in place.
report_dir <- if (file.exists(target_report) || !file.exists(legacy_report)) full_dir else legacy_report_dir
fig_dir <- file.path(report_dir, "figures")

paths <- list(
  report_html = if (identical(report_dir, full_dir)) target_report else legacy_report,
  evidence_summary = file.path(report_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fibrtraq_evidence-summary.md"),
  result_csv = file.path(full_dir, "sub-dtor985_desc-dtor985-left-thalamus-harvardoxford-exploratory-full_fibrtraq_per-voxel-hits.csv"),
  rsc_csv = file.path(full_dir, "sub-dtor985_desc-dtor985-left-thalamus-harvardoxford-exploratory-full_fibrtraq_rsc-summary.csv"),
  validation_json = file.path(full_dir, "sub-dtor985_desc-dtor985-left-thalamus-harvardoxford-exploratory-full_fibrtraq_validation.json"),
  exploratory_csv = file.path(report_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fibrtraq_exploratory-per-voxel.csv"),
  target_totals_csv = file.path(report_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fibrtraq_exploratory-target-totals.csv"),
  sampled_streamline_points_csv = file.path(report_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_sampled-streamline-points.csv"),
  sampled_streamline_summary_csv = file.path(report_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_sampled-streamline-summary.csv"),
  figures = list(
    target_totals = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_exploratory-target-totals.png"),
    primary_target = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_exploratory-primary-target.png"),
    target_diversity = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_exploratory-target-diversity.png"),
    target_entropy = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_exploratory-target-entropy.png"),
    leave_no_target = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_exploratory-leave-no-target-fraction.png"),
    sampled_streamlines = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines.png"),
    sampled_streamlines_3d = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines-3d-interactive.html"),
    sampled_streamlines_3d_snapshot = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines-3d.png"),
    sampled_streamlines_anatomical = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines-anatomical-overlay.png")
  )
)

# Prefer the landing page; fall back to a generated report or viewer.
entrypoints <- c(
  file.path(left_root, "fibrtraq-report-index-portable.html"),
  file.path(left_root, "fibrtraq-report-index.html"),
  paths$report_html,
  paths$figures$sampled_streamlines_3d,
  legacy_report,
  file.path(legacy_report_dir, "figures", paste0(report_prefix, "_fig_sampled-streamlines-3d-interactive.html"))
)
available_entrypoints <- entrypoints[file.exists(entrypoints)]
if (!length(available_entrypoints)) {
  stop("No results index, report or 3-D viewer was found under: ", left_root,
       ". Check FIBRTRAQ_DTOR_ROOT and extract the complete results folder.",
       call. = FALSE)
}
open_output <- function(path) {
  FibrTraQ:::ftq_qc_open_output(path)
}
cat("\nOpening results:\n", available_entrypoints[1L], "\n", sep = "")
open_output(available_entrypoints[1L])

# These optional console summaries do not prevent access to completed reports.
summary_inputs <- unlist(paths[c("result_csv", "validation_json")], use.names = FALSE)
if (all(file.exists(summary_inputs))) tryCatch({
  result <- read.csv(paths$result_csv, check.names = FALSE)
  if (all(file.exists(c(paths$exploratory_csv, paths$target_totals_csv)))) {
    exploratory <- read.csv(paths$exploratory_csv, check.names = FALSE)
    target_totals <- read.csv(paths$target_totals_csv, check.names = FALSE)
  } else {
    # Derive the optional console summary in memory; no extra CSVs or maps.
    derived <- FibrTraQ::ftq_exploratory_summary(paths$result_csv, write_maps = FALSE)
    exploratory <- derived$data
    target_totals <- derived$target_totals
  }
  validation <- jsonlite::fromJSON(paths$validation_json, simplifyVector = FALSE)
  target_association_col <- if ("focal_voxel_streamline_target_associations" %in% names(target_totals)) {
    "focal_voxel_streamline_target_associations"
  } else {
    "tracks"
  }

  connected <- result$tracks > 0
  summary_table <- data.frame(
    metric = c(
      "validation_passed",
      "focal_voxels",
      "focal_voxels_with_tracks",
      "focal_voxel_streamline_associations",
      "leave_no_target_associations",
      "focal_voxel_streamline_target_associations",
      "nonzero_targets",
      "median_target_diversity_connected_voxels",
      "median_entropy_fraction_connected_voxels"
    ),
    value = c(
      as.character(validation$passed),
      format(nrow(result), big.mark = ","),
      paste0(format(sum(connected), big.mark = ","), " / ", format(nrow(result), big.mark = ",")),
      format(sum(result$tracks), big.mark = ","),
      format(sum(result$leave_no_target), big.mark = ","),
      format(sum(target_totals[[target_association_col]]), big.mark = ","),
      paste0(format(sum(target_totals[[target_association_col]] > 0), big.mark = ","), " / ", format(nrow(target_totals), big.mark = ",")),
      format(stats::median(exploratory$target_diversity[exploratory$tracks > 0]), big.mark = ","),
      sprintf("%.3f", stats::median(exploratory$target_entropy_fraction[exploratory$tracks > 0]))
    ),
    stringsAsFactors = FALSE
  )

  cat("\nFibrTraQ dTOR-985 full validation summary\n")
  print(summary_table, row.names = FALSE)

  cat("\nTop 20 target association totals\n")
  print(utils::head(target_totals, 20), row.names = FALSE)

}, error = function(e) {
  message("The report opened, but the optional console summary is unavailable: ",
          conditionMessage(e))
}) else {
  cat("Supporting tables are incomplete; the available report was opened.\n")
}

cat("\nFigure paths:\n")
print(paths$figures[file.exists(unlist(paths$figures, use.names = FALSE))])

cat("\nOpen a specific figure from RStudio with, for example:\n")
cat("open_output(paths$figures$primary_target)\n")
cat("open_output(paths$figures$sampled_streamlines_3d)\n")
cat("open_output(paths$figures$sampled_streamlines_anatomical)\n")
cat("open_output(paths$figures$sampled_streamlines)\n")
cat("open_output(paths$figures$target_diversity)\n")
cat("open_output(paths$figures$target_entropy)\n")
