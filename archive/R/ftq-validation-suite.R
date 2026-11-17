#' Run a suite of FibrTraQ validation manifests
#'
#' Run preflight checks, and optionally full FibrTraQ analyses, across one or
#' more manifest files. This is intended for publication validation where each
#' dataset/application is represented by a manifest and the suite output becomes
#' an auditable summary table.
#'
#' @param manifests Character vector of manifest paths, or a directory containing
#'   `.yml`, `.yaml`, or `.json` manifests. Directories are searched
#'   recursively.
#' @param output_csv Optional path where the aggregate suite table should be
#'   written.
#' @param preflight_only If TRUE, only run manifest preflight checks.
#' @param overwrite Whether existing FibrTraQ outputs may be overwritten during
#'   full runs.
#' @param write_maps Optional override for manifest map output settings.
#' @param write_report Optional override for manifest report output settings.
#' @param continue_on_error If TRUE, failed manifests are retained as rows in the
#'   output table. If FALSE, the first failed manifest stops the suite.
#' @param verbose Emit progress messages from the underlying FibrTraQ run.
#'
#' @return A data.frame with one row per manifest and summary columns for
#'   preflight status, validation status, core output paths, and run totals.
#' @export
ftq_run_validation_suite <- function(manifests,
                                     output_csv = NULL,
                                     preflight_only = FALSE,
                                     overwrite = FALSE,
                                     write_maps = NULL,
                                     write_report = NULL,
                                     continue_on_error = TRUE,
                                     verbose = FALSE) {
  manifest_paths <- ftq_validation_manifest_paths(manifests)
  rows <- lapply(manifest_paths, function(path) {
    ftq_run_validation_manifest(
      path = path,
      preflight_only = preflight_only,
      overwrite = overwrite,
      write_maps = write_maps,
      write_report = write_report,
      continue_on_error = continue_on_error,
      verbose = verbose
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  if (!is.null(output_csv)) {
    output_csv <- path.expand(output_csv)
    dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, output_csv, row.names = FALSE)
  }

  out
}

ftq_validation_manifest_paths <- function(manifests) {
  if (is.null(manifests) || !length(manifests)) {
    stop("manifests must be a manifest path, a vector of paths, or a directory.")
  }

  if (is.character(manifests) && length(manifests) == 1L &&
      dir.exists(path.expand(manifests))) {
    root <- path.expand(manifests)
    paths <- list.files(
      root,
      pattern = "\\.(json|ya?ml)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  } else {
    paths <- path.expand(as.character(manifests))
  }

  paths <- paths[nzchar(paths)]
  if (!length(paths)) stop("No manifest files were found.")
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Manifest file(s) do not exist:\n", paste(missing, collapse = "\n"))
  }

  normalizePath(paths, mustWork = TRUE)
}

ftq_run_validation_manifest <- function(path,
                                        preflight_only,
                                        overwrite,
                                        write_maps,
                                        write_report,
                                        continue_on_error,
                                        verbose) {
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  m <- tryCatch(ftq_read_manifest(path), error = function(e) e)
  if (inherits(m, "error")) {
    if (!isTRUE(continue_on_error)) stop(m)
    return(ftq_validation_suite_row(
      manifest_path = path,
      run_status = "manifest_error",
      started_at = started,
      completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      error = conditionMessage(m)
    ))
  }

  preflight_args <- list(manifest = m, overwrite = overwrite)
  if (!is.null(write_maps)) preflight_args$write_maps <- isTRUE(write_maps)
  preflight <- tryCatch(do.call(ftq_preflight, preflight_args), error = function(e) e)
  if (inherits(preflight, "error")) {
    if (!isTRUE(continue_on_error)) stop(preflight)
    return(ftq_validation_suite_row(
      manifest = m,
      manifest_path = path,
      run_status = "preflight_error",
      started_at = started,
      completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      error = conditionMessage(preflight)
    ))
  }

  base_row <- ftq_validation_suite_row(
    manifest = m,
    manifest_path = path,
    preflight = preflight,
    run_status = if (isTRUE(preflight_only)) "preflight_only" else "preflight_passed",
    started_at = started
  )

  if (isTRUE(preflight_only)) {
    base_row$completed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    return(base_row)
  }

  if (!isTRUE(preflight$passed)) {
    if (!isTRUE(continue_on_error)) {
      stop("FibrTraQ preflight failed for ", path, ":\n",
           paste(ftq_preflight_failures(preflight), collapse = "\n"))
    }
    base_row$run_status <- "preflight_failed"
    base_row$completed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    return(base_row)
  }

  run_args <- list(
    manifest = m,
    overwrite = overwrite,
    verbose = verbose
  )
  if (!is.null(write_maps)) run_args$write_maps <- isTRUE(write_maps)
  if (!is.null(write_report)) run_args$write_report <- isTRUE(write_report)

  run <- tryCatch(do.call(ftq_run, run_args), error = function(e) e)
  if (inherits(run, "error")) {
    if (!isTRUE(continue_on_error)) stop(run)
    base_row$run_status <- "run_error"
    base_row$error <- conditionMessage(run)
    base_row$completed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    return(base_row)
  }

  ftq_validation_suite_row(
    manifest = run$manifest,
    manifest_path = path,
    preflight = run$preflight,
    validation = run$validation,
    outputs = run$outputs,
    report = run$report,
    result = run$result,
    rsc = run$rsc,
    run_status = "completed",
    started_at = started,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
}

ftq_validation_suite_row <- function(manifest = NULL,
                                     manifest_path = NA_character_,
                                     preflight = NULL,
                                     validation = NULL,
                                     outputs = NULL,
                                     report = NULL,
                                     result = NULL,
                                     rsc = NULL,
                                     run_status,
                                     started_at,
                                     completed_at = NA_character_,
                                     error = "") {
  preflight_counts <- ftq_preflight_counts(preflight)
  validation_counts <- ftq_validation_counts(validation)
  summary <- ftq_suite_result_summary(result, rsc)

  data.frame(
    manifest_path = as.character(manifest_path),
    dataset_root = ftq_suite_field(manifest$dataset_root),
    subject = ftq_suite_field(manifest$subject),
    session = ftq_suite_field(manifest$session),
    description = ftq_suite_field(manifest$description),
    prefix = ftq_suite_field(manifest$prefix),
    output_dir = ftq_suite_field(manifest$output_dir),
    run_status = run_status,
    preflight_passed = ftq_suite_logical(preflight$passed),
    preflight_errors = preflight_counts$errors,
    preflight_warnings = preflight_counts$warnings,
    preflight_recommendations = preflight_counts$recommendations,
    validation_passed = ftq_suite_logical(validation$passed),
    validation_failed_checks = validation_counts$failed_checks,
    endpoint_analysis = ftq_suite_field(manifest$endpoint_analysis),
    biological_qc_status = ftq_suite_field(validation$biological_qc_status),
    focal_voxels = summary$focal_voxels,
    focal_voxel_streamline_associations = summary$focal_voxel_streamline_associations,
    local_associations = summary$local_associations,
    leave_no_target_associations = summary$leave_no_target_associations,
    leave_with_target_associations = summary$leave_with_target_associations,
    focal_voxel_streamline_target_associations = summary$focal_voxel_streamline_target_associations,
    result_csv = ftq_suite_field(outputs$result_csv),
    rsc_csv = ftq_suite_field(outputs$rsc_csv),
    fate_csv = ftq_suite_field(outputs$fate_csv),
    roi_summary_csv = ftq_suite_field(outputs$roi_summary_csv),
    endpoint_summary_csv = ftq_suite_field(outputs$endpoint_summary_csv),
    endpoint_penetration_csv = ftq_suite_field(outputs$endpoint_penetration_csv),
    endpoint_context_csv = ftq_suite_field(outputs$endpoint_context_csv),
    metadata_json = ftq_suite_field(outputs$metadata_json),
    validation_json = ftq_suite_field(outputs$validation_json),
    report_html = ftq_suite_field(report$html),
    error = as.character(error),
    started_at = started_at,
    completed_at = completed_at,
    target_intersection_mode =
      ftq_suite_field(manifest$target_intersection_mode),
    local_search_radius_mm = suppressWarnings(as.numeric(
      ftq_suite_field(manifest$local_search_radius_mm)
    )),
    stringsAsFactors = FALSE
  )
}

ftq_suite_field <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return(NA_character_)
  as.character(x[1])
}

ftq_suite_logical <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return(NA)
  isTRUE(x[1])
}

ftq_preflight_counts <- function(preflight) {
  if (is.null(preflight) || is.null(preflight$checks) || !nrow(preflight$checks)) {
    return(list(errors = NA_integer_, warnings = NA_integer_, recommendations = NA_integer_))
  }
  failed <- !preflight$checks$passed
  list(
    errors = sum(failed & preflight$checks$severity == "error"),
    warnings = sum(failed & preflight$checks$severity == "warning"),
    recommendations = sum(failed & preflight$checks$severity == "recommendation")
  )
}

ftq_validation_counts <- function(validation) {
  if (is.null(validation) || is.null(validation$checks) || !nrow(validation$checks)) {
    return(list(failed_checks = NA_character_))
  }
  failed <- validation$checks$check[!validation$checks$passed]
  list(failed_checks = paste(failed, collapse = ";"))
}

ftq_suite_result_summary <- function(result, rsc) {
  if (is.null(result) || is.null(rsc)) {
    return(list(
      focal_voxels = NA_integer_,
      focal_voxel_streamline_associations = NA_real_,
      local_associations = NA_real_,
      leave_no_target_associations = NA_real_,
      leave_with_target_associations = NA_real_,
      focal_voxel_streamline_target_associations = NA_character_
    ))
  }

  summary <- ftq_result_summary(result, rsc)
  target_associations <- summary$voxelwise_target_associations
  target_detail <- if (length(target_associations)) {
    paste(paste0(names(target_associations), "=", unlist(target_associations, use.names = FALSE)), collapse = ";")
  } else {
    ""
  }

  list(
    focal_voxels = summary$focal_voxels,
    focal_voxel_streamline_associations = summary$voxelwise_track_associations,
    local_associations = summary$voxelwise_local_associations,
    leave_no_target_associations = summary$voxelwise_leave_no_target_associations,
    leave_with_target_associations = summary$voxelwise_leave_with_target_associations,
    focal_voxel_streamline_target_associations = target_detail
  )
}
