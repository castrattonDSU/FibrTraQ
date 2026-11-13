#' Locate the bundled FibrTraQ known-answer example
#'
#' Return absolute paths to the package-owned synthetic tractogram, masks,
#' reference image, and seven expected scientific output tables used by
#' [ftq_run_example()]. The fixture is intentionally small and is not
#' anatomical MRI or evidence of biological validity.
#'
#' By default, the function verifies the size and MD5 checksum of every bundled
#' input and expected table before returning paths. The returned claim registry
#' defines 18 independently specified known answers and the exact table fields
#' from which they must be recovered.
#'
#' @param verify Whether to verify the bundled file inventory and checksums.
#'
#' @return A list containing the fixture root, input paths, expected-table
#'   paths, analysis parameters, known-answer claims, and checksum inventory.
#' @export
#'
#' @examples
#' example <- ftq_example_data()
#' example$tractogram
#' example$target_masks
ftq_example_data <- function(verify = TRUE) {
  if (!is.logical(verify) || length(verify) != 1L || is.na(verify)) {
    stop("verify must be TRUE or FALSE.", call. = FALSE)
  }
  root <- system.file(
    "extdata",
    "fibrtraq-known-answer",
    package = "FibrTraQ"
  )
  if (!nzchar(root) || !dir.exists(root)) {
    stop(
      "The bundled FibrTraQ known-answer example is unavailable. ",
      "Reinstall FibrTraQ from a complete package build.",
      call. = FALSE
    )
  }
  root <- normalizePath(root, mustWork = TRUE)

  checksum_path <- file.path(root, "checksums.csv")
  inventory <- utils::read.csv(
    checksum_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_inventory <- c("relative_path", "size_bytes", "md5")
  if (!all(required_inventory %in% names(inventory)) ||
      !nrow(inventory) || anyDuplicated(inventory$relative_path)) {
    stop("The bundled example checksum inventory is malformed.",
         call. = FALSE)
  }
  expected_inventory <- c(
    "analysis-parameters.json",
    "known-answer-claims.csv",
    "README.md",
    file.path(
      "inputs",
      c(
        "sub-example_context.nii.gz",
        "sub-example_focal.nii.gz",
        "sub-example_local.nii.gz",
        "sub-example_reference.nii.gz",
        "sub-example_target-brainstem.nii.gz",
        "sub-example_target-cortex.nii.gz",
        "sub-example_tractogram.trk"
      )
    ),
    file.path(
      "expected",
      paste0(
        c(
          "endpoint_context",
          "endpoint_penetration",
          "endpoint_summary",
          "fate",
          "result",
          "roi_summary",
          "rsc"
        ),
        ".csv"
      )
    )
  )
  if (!setequal(inventory$relative_path, expected_inventory)) {
    missing <- setdiff(expected_inventory, inventory$relative_path)
    extra <- setdiff(inventory$relative_path, expected_inventory)
    detail <- c(
      if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(extra)) paste0("unexpected: ", paste(extra, collapse = ", "))
    )
    stop(
      "The bundled example checksum inventory is incomplete: ",
      paste(detail, collapse = "; "),
      ".",
      call. = FALSE
    )
  }
  inventory$path <- file.path(root, inventory$relative_path)

  if (isTRUE(verify)) {
    exists <- file.exists(inventory$path)
    observed_size <- rep(NA_real_, nrow(inventory))
    observed_md5 <- rep(NA_character_, nrow(inventory))
    observed_size[exists] <- as.numeric(file.info(inventory$path[exists])$size)
    observed_md5[exists] <- unname(tools::md5sum(inventory$path[exists]))
    valid <- exists &
      observed_size == inventory$size_bytes &
      observed_md5 == inventory$md5
    if (!all(valid)) {
      bad <- inventory$relative_path[!valid]
      stop(
        "The bundled FibrTraQ example failed its integrity check: ",
        paste(bad, collapse = ", "),
        ". Reinstall the package before using the example.",
        call. = FALSE
      )
    }
  }

  input <- function(name) {
    path <- file.path(root, "inputs", name)
    if (!file.exists(path)) {
      stop("Bundled example input is missing: ", name, call. = FALSE)
    }
    normalizePath(path, mustWork = TRUE)
  }
  expected <- function(name) {
    path <- file.path(root, "expected", paste0(name, ".csv"))
    if (!file.exists(path)) {
      stop("Bundled expected table is missing: ", name, call. = FALSE)
    }
    normalizePath(path, mustWork = TRUE)
  }

  parameters <- jsonlite::read_json(
    file.path(root, "analysis-parameters.json"),
    simplifyVector = TRUE
  )
  required_parameters <- c(
    "dataset_id",
    "dataset_type",
    "biological_interpretation",
    "subject",
    "focal_label",
    "target_labels",
    "description",
    "track_space",
    "discretization",
    "axis_flip",
    "focal_assignment_mode",
    "streamline_voxelization_mode",
    "target_intersection_mode",
    "endpoint_analysis",
    "endpoint_context_labels",
    "track_reader",
    "track_chunk_size",
    "expected_table_backend",
    "example_execution",
    "expected_table_generation",
    "generation_provenance"
  )
  if (!all(required_parameters %in% names(parameters)) ||
      !identical(parameters$dataset_id, "fibrtraq-known-answer-v1") ||
      !identical(parameters$biological_interpretation, FALSE) ||
      !identical(parameters$target_intersection_mode, "full_streamline") ||
      !identical(parameters$expected_table_backend, "r")) {
    stop("The bundled example parameter record is malformed.",
         call. = FALSE)
  }
  valid_logical_fields <- function(x, fields) {
    is.list(x) &&
      all(fields %in% names(x)) &&
      all(vapply(x[fields], function(value) {
        is.logical(value) && length(value) == 1L && !is.na(value)
      }, logical(1)))
  }
  example_flags <- c(
    "write_maps",
    "write_report",
    "write_exploratory_report",
    "write_qc_viewer",
    "progress"
  )
  generation_flags <- c(
    "write_maps",
    "write_report",
    "write_exploratory_report",
    "write_qc_viewer",
    "progress"
  )
  provenance_fields <- c(
    "generator",
    "package_version",
    "r_version",
    "RNifti_version",
    "jsonlite_version",
    "devtools_version"
  )
  provenance <- parameters$generation_provenance
  example_flag_values <- if (is.list(parameters$example_execution)) {
    vapply(
      parameters$example_execution[example_flags],
      isTRUE,
      logical(1)
    )
  } else {
    logical()
  }
  generation_flag_values <-
    if (is.list(parameters$expected_table_generation)) {
      vapply(
        parameters$expected_table_generation[generation_flags],
        isTRUE,
        logical(1)
      )
    } else {
      logical()
    }
  if (!valid_logical_fields(parameters$example_execution, example_flags) ||
      !valid_logical_fields(
        parameters$expected_table_generation,
        generation_flags
      ) ||
      !identical(
        unname(example_flag_values),
        c(TRUE, TRUE, TRUE, TRUE, FALSE)
      ) ||
      !identical(
        unname(generation_flag_values),
        rep(FALSE, length(generation_flags))
      ) ||
      !identical(
        as.integer(parameters$example_execution$sample_per_target),
        6L
      ) ||
      !is.list(provenance) ||
      !all(provenance_fields %in% names(provenance)) ||
      any(!nzchar(vapply(
        provenance[provenance_fields],
        as.character,
        character(1)
      )))) {
    stop("The bundled example parameter record is malformed.",
         call. = FALSE)
  }
  claims <- utils::read.csv(
    file.path(root, "known-answer-claims.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  claim_columns <- c(
    "claim",
    "source_table",
    "value_column",
    "filter_column_1",
    "filter_value_1",
    "filter_column_2",
    "filter_value_2",
    "expected_value"
  )
  if (!identical(names(claims), claim_columns) ||
      !nrow(claims) ||
      anyDuplicated(claims$claim)) {
    stop("The bundled example known-answer claims are malformed.",
         call. = FALSE)
  }

  missing_text <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))
  valid_claim_sources <- c(
    "roi_summary",
    "endpoint_summary",
    "rsc",
    "endpoint_context",
    "fate"
  )
  filter_1_incomplete <- xor(
    missing_text(claims$filter_column_1),
    missing_text(claims$filter_value_1)
  )
  filter_2_incomplete <- xor(
    missing_text(claims$filter_column_2),
    missing_text(claims$filter_value_2)
  )
  if (any(missing_text(claims$claim)) ||
      any(missing_text(claims$source_table)) ||
      any(missing_text(claims$value_column)) ||
      any(!claims$source_table %in% valid_claim_sources) ||
      any(filter_1_incomplete) ||
      any(filter_2_incomplete) ||
      any(!is.finite(claims$expected_value))) {
    stop("The bundled example known-answer claims are malformed.",
         call. = FALSE)
  }

  list(
    id = "fibrtraq-known-answer-v1",
    root = root,
    tractogram = input("sub-example_tractogram.trk"),
    focal_mask = input("sub-example_focal.nii.gz"),
    target_masks = c(
      cortex = input("sub-example_target-cortex.nii.gz"),
      brainstem = input("sub-example_target-brainstem.nii.gz")
    ),
    local_mask = input("sub-example_local.nii.gz"),
    reference_image = input("sub-example_reference.nii.gz"),
    endpoint_contexts = c(
      analysis_domain = input("sub-example_context.nii.gz")
    ),
    expected = c(
      result = expected("result"),
      rsc = expected("rsc"),
      fate = expected("fate"),
      roi_summary = expected("roi_summary"),
      endpoint_summary = expected("endpoint_summary"),
      endpoint_penetration = expected("endpoint_penetration"),
      endpoint_context = expected("endpoint_context")
    ),
    parameters = parameters,
    claims = claims,
    inventory = inventory
  )
}

#' Run and verify the bundled FibrTraQ example
#'
#' Run the complete simple file workflow on the bundled known-answer
#' tractogram and masks. The function writes the standard and exploratory
#' reports, NIfTI maps, endpoint tables, and interactive 3-D viewer, then
#' compares seven generated scientific CSV tables with the versioned expected
#' tables installed with the package. It also resolves all 18 independently
#' declared known-answer claims from both the stored and newly generated tables.
#'
#' The comparison is exact after normalizing two provenance-only fields:
#' backend identity and the installation-specific absolute path to the endpoint
#' context mask. All scientific values, labels, denominators, columns, rows,
#' and ordering must match. Four HTML files are also checked for existence and
#' nonzero size. These generation checks do not verify browser rendering or
#' biological validity. The results index displays both types of evidence.
#'
#' @param output_dir Directory for generated example outputs. The default is a
#'   new temporary directory.
#' @param backend Counting backend passed to [ftq_run_files()].
#' @param overwrite Whether an existing `output_dir` may be replaced.
#' @param verbose Emit run messages.
#' @param open Open the report index after successful verification.
#'
#' @return An object of class `ftq_example_result` containing the bundled data
#'   paths, the [ftq_run_files()] result, exact table and claim comparison
#'   details, and report paths.
#' @export
#'
#' @examples
#' \donttest{
#' example <- ftq_run_example(verbose = FALSE)
#' example$validation$passed
#' example$report_index
#' }
ftq_run_example <- function(
    output_dir = tempfile("ftq-"),
    backend = c("auto", "cpp", "r"),
    overwrite = FALSE,
    verbose = TRUE,
    open = FALSE) {
  backend <- match.arg(backend)
  if (length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(trimws(output_dir))) {
    stop("output_dir must be one non-empty path.", call. = FALSE)
  }
  output_dir <- normalizePath(
    path.expand(output_dir),
    mustWork = FALSE
  )

  example <- ftq_example_data(verify = TRUE)
  params <- example$parameters
  analysis <- ftq_run_files(
    tractogram = example$tractogram,
    focal_mask = example$focal_mask,
    target_masks = example$target_masks,
    output_dir = output_dir,
    subject = params$subject,
    focal_label = params$focal_label,
    local_mask = example$local_mask,
    reference_image = example$reference_image,
    endpoint_contexts = example$endpoint_contexts,
    description = params$description,
    track_space = params$track_space,
    discretization = params$discretization,
    axis_flip = as.logical(params$axis_flip),
    focal_assignment_mode = params$focal_assignment_mode,
    streamline_voxelization_mode = params$streamline_voxelization_mode,
    assignment_mode = params$target_intersection_mode,
    endpoint_analysis = params$endpoint_analysis,
    track_reader = params$track_reader,
    track_chunk_size = as.integer(params$track_chunk_size),
    backend = backend,
    write_maps = isTRUE(params$example_execution$write_maps),
    write_report = isTRUE(params$example_execution$write_report),
    write_exploratory_report = isTRUE(
      params$example_execution$write_exploratory_report
    ),
    write_qc_viewer = isTRUE(params$example_execution$write_qc_viewer),
    sample_per_target = as.integer(
      params$example_execution$sample_per_target
    ),
    overwrite = overwrite,
    verbose = verbose,
    progress = isTRUE(params$example_execution$progress),
    open = FALSE
  )

  run <- analysis$index[1, , drop = FALSE]
  actual <- c(
    result = as.character(run$result_csv),
    rsc = as.character(run$rsc_csv),
    fate = as.character(run$fate_csv),
    roi_summary = as.character(run$roi_summary_csv),
    endpoint_summary = as.character(run$endpoint_summary_csv),
    endpoint_penetration = as.character(run$endpoint_penetration_csv),
    endpoint_context = as.character(run$endpoint_context_csv)
  )
  tables <- ftq_example_compare_tables(actual, example$expected)
  claims <- ftq_example_check_claims(
    actual = actual,
    expected = example$expected,
    claims = example$claims
  )

  report_paths <- c(
    report_index = analysis$report_index,
    primary_report = as.character(run$primary_report),
    exploratory_report = as.character(run$exploratory_report),
    interactive_3d_qc = as.character(run$interactive_3d_qc)
  )
  report_size <- as.numeric(file.info(report_paths)$size)
  report_ok <- file.exists(report_paths) &
    !is.na(report_size) &
    report_size > 0
  reports <- data.frame(
    artifact = names(report_paths),
    path = unname(report_paths),
    passed = unname(report_ok),
    details = ifelse(report_ok, "generated and non-empty", "missing or empty"),
    stringsAsFactors = FALSE
  )
  passed <- isTRUE(run$validation_passed) &&
    all(tables$passed) &&
    all(claims$passed) &&
    all(reports$passed)
  validation <- list(
    passed = passed,
    package_validation_passed = isTRUE(run$validation_passed),
    tables = tables,
    claims = claims,
    reports = reports
  )
  # Add observed verification results only to this invocation's new index.
  ftq_write_regional_study_index(
    run_index = analysis$index,
    path = analysis$report_index,
    output_dir = output_dir,
    description = params$description,
    example_validation = validation
  )

  if (!passed) {
    failed <- c(
      tables$artifact[!tables$passed],
      paste0("claim:", claims$claim[!claims$passed]),
      reports$artifact[!reports$passed],
      if (!isTRUE(run$validation_passed)) "package_validation"
    )
    stop(
      "The bundled FibrTraQ example failed verification: ",
      paste(failed, collapse = ", "),
      ". Outputs were preserved at ", output_dir, ".",
      call. = FALSE
    )
  }

  out <- structure(
    list(
      data = example,
      analysis = analysis,
      validation = validation,
      output_dir = output_dir,
      report_index = analysis$report_index,
      primary_report = as.character(run$primary_report),
      exploratory_report = as.character(run$exploratory_report),
      interactive_3d_qc = as.character(run$interactive_3d_qc)
    ),
    class = "ftq_example_result"
  )
  if (isTRUE(open)) ftq_qc_open_output(out$report_index)
  out
}

#' @export
print.ftq_example_result <- function(x, ...) {
  passed_tables <- sum(x$validation$tables$passed)
  total_tables <- nrow(x$validation$tables)
  passed_claims <- sum(x$validation$claims$passed)
  total_claims <- nrow(x$validation$claims)
  passed_reports <- sum(x$validation$reports$passed)
  total_reports <- nrow(x$validation$reports)
  cat(
    "FibrTraQ known-answer example (synthetic data)\n",
    "  Exact tables: ", passed_tables, "/", total_tables, "\n",
    "  Known answers: ", passed_claims, "/", total_claims, "\n",
    "  HTML files:   ", passed_reports, "/", total_reports,
    " generated and non-empty (browser rendering not checked)\n",
    "  Computational validation: ",
    if (isTRUE(x$validation$package_validation_passed)) "passed" else "failed", "\n",
    "  Example verification: ",
    if (isTRUE(x$validation$passed)) "passed" else "failed", "\n",
    "  Biological validity is not assessed by this synthetic example.\n",
    "  Output:       ", x$output_dir, "\n",
    "  Report index: ", x$report_index, "\n",
    sep = ""
  )
  invisible(x)
}

ftq_example_compare_tables <- function(actual, expected) {
  artifacts <- names(expected)
  if (is.null(artifacts) || !length(artifacts) ||
      !identical(names(actual), artifacts)) {
    stop("Known-answer table mappings are inconsistent.", call. = FALSE)
  }

  rows <- lapply(artifacts, function(artifact) {
    actual_path <- actual[[artifact]]
    expected_path <- expected[[artifact]]
    if (!file.exists(actual_path) || !file.exists(expected_path)) {
      return(data.frame(
        artifact = artifact,
        actual_rows = NA_integer_,
        expected_rows = NA_integer_,
        passed = FALSE,
        details = "actual or expected table is missing",
        stringsAsFactors = FALSE
      ))
    }
    observed <- ftq_example_read_table(actual_path, artifact)
    known <- ftq_example_read_table(expected_path, artifact)
    passed <- identical(observed, known)
    detail <- if (passed) {
      "exact scientific-table match"
    } else {
      paste(
        utils::head(
          as.character(all.equal(observed, known, tolerance = 0)),
          3L
        ),
        collapse = "; "
      )
    }
    data.frame(
      artifact = artifact,
      actual_rows = nrow(observed),
      expected_rows = nrow(known),
      passed = passed,
      details = detail,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ftq_example_check_claims <- function(actual, expected, claims) {
  if (!identical(names(actual), names(expected))) {
    stop("Known-answer table mappings are inconsistent.", call. = FALSE)
  }
  read_tables <- function(paths) {
    stats::setNames(
      lapply(names(paths), function(artifact) {
        path <- paths[[artifact]]
        if (!file.exists(path)) return(NULL)
        ftq_example_read_table(path, artifact)
      }),
      names(paths)
    )
  }
  actual_tables <- read_tables(actual)
  expected_tables <- read_tables(expected)

  extract_value <- function(tables, claim) {
    source <- as.character(claim$source_table)
    table <- tables[[source]]
    if (is.null(table)) {
      return(list(value = NA_real_, error = "source table is missing"))
    }
    value_column <- as.character(claim$value_column)
    if (!value_column %in% names(table)) {
      return(list(
        value = NA_real_,
        error = paste0("value column is missing: ", value_column)
      ))
    }
    keep <- rep(TRUE, nrow(table))
    for (suffix in c("1", "2")) {
      filter_column <- claim[[paste0("filter_column_", suffix)]]
      filter_value <- claim[[paste0("filter_value_", suffix)]]
      has_filter <- !is.na(filter_column) && nzchar(trimws(filter_column))
      if (!has_filter) next
      if (!filter_column %in% names(table)) {
        return(list(
          value = NA_real_,
          error = paste0("filter column is missing: ", filter_column)
        ))
      }
      keep <- keep & as.character(table[[filter_column]]) == filter_value
    }
    selected <- which(keep)
    if (length(selected) != 1L) {
      return(list(
        value = NA_real_,
        error = paste0(
          "claim selector matched ",
          length(selected),
          " rows"
        )
      ))
    }
    value <- suppressWarnings(as.numeric(table[[value_column]][selected]))
    if (length(value) != 1L || !is.finite(value)) {
      return(list(value = NA_real_, error = "claim value is not finite numeric"))
    }
    list(value = value, error = NA_character_)
  }

  rows <- lapply(seq_len(nrow(claims)), function(i) {
    claim <- claims[i, , drop = FALSE]
    stored <- extract_value(expected_tables, claim)
    observed <- extract_value(actual_tables, claim)
    declared <- as.numeric(claim$expected_value)
    stored_passed <- is.na(stored$error) &&
      isTRUE(all.equal(stored$value, declared, tolerance = 0))
    actual_passed <- is.na(observed$error) &&
      isTRUE(all.equal(observed$value, declared, tolerance = 0))
    passed <- stored_passed && actual_passed
    detail <- if (passed) {
      "stored and generated values match the declared known answer"
    } else {
      paste(
        c(
          if (!stored_passed) {
            paste0(
              "stored=",
              if (is.na(stored$error)) stored$value else stored$error
            )
          },
          if (!actual_passed) {
            paste0(
              "generated=",
              if (is.na(observed$error)) observed$value else observed$error
            )
          },
          paste0("declared=", declared)
        ),
        collapse = "; "
      )
    }
    data.frame(
      claim = as.character(claim$claim),
      source_table = as.character(claim$source_table),
      expected_value = declared,
      stored_value = stored$value,
      actual_value = observed$value,
      stored_passed = stored_passed,
      actual_passed = actual_passed,
      passed = passed,
      details = detail,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ftq_example_read_table <- function(path, artifact) {
  table <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (identical(artifact, "roi_summary") && "backend" %in% names(table)) {
    table$backend <- "<normalized-backend-provenance>"
  }
  if (identical(artifact, "endpoint_context") &&
      "source_path" %in% names(table)) {
    has_source <- !is.na(table$source_path) & nzchar(table$source_path)
    table$source_path[has_source] <- basename(table$source_path[has_source])
  }
  table
}
