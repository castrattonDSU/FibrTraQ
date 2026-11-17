# Check the same prepared inputs and per-run manifests used by the real analysis.
# Scratch files are needed only for inputs such as a requested union or a 4-D reference.
ftq_check_study_inputs <- function(inputs, focal_regions, output_dir, check_destination,
                                   prepare_subject, make_manifest, overwrite, check_run_destination) {
  scratch <- tempfile("ftq-check-")
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  failure <- function(check, detail) {
    data.frame(check = check, severity = "error", passed = FALSE, detail = detail,
               stringsAsFactors = FALSE)
  }
  output_error <- function(e) failure("output_directory", conditionMessage(e))
  destination_error <- tryCatch(check_destination(), error = output_error)

  rows <- list()
  for (subject in inputs$subjects) {
    masks <- inputs$region_masks[[subject]]
    prep_paths <- c(stats::setNames(unname(masks), paste0("mask:", names(masks))),
                    reference_image = inputs$reference_images[[subject]],
                    background_image = inputs$background_images[[subject]],
                    local_mask = inputs$local_masks[[subject]])
    if (is.character(inputs$bvals[[subject]])) {
      prep_paths <- c(prep_paths, bvals = inputs$bvals[[subject]])
    }
    trk <- inputs$tractograms[[subject]]
    contexts <- unlist(inputs$endpoint_contexts[[subject]], use.names = TRUE)
    if (length(contexts)) names(contexts) <- paste0("endpoint_context:", names(contexts))
    paths <- c(prep_paths, tractogram = unname(trk), contexts)
    bad <- stats::setNames(!file.exists(paths) | dir.exists(paths) | file.access(paths, 4L) != 0L,
                           names(paths))
    errors <- if (any(bad)) {
      failure(paste0("readable_file:", names(paths)[bad]), paste("Cannot read", paths[bad]))
    }
    # Open the TRK header through the same reader and signature checks as counting.
    # TCK headers are checked by ftq_preflight() for each coordinate reference.
    if (!bad[["tractogram"]] && grepl("[.]trk$", trk, ignore.case = TRUE)) {
      header_error <- tryCatch({ftq_read_tractogram(trk, reader = "stream"); NULL},
                               error = function(e) failure("trk_header_readable",
                                                           paste(trk, conditionMessage(e))))
      errors <- rbind(errors, header_error)
    }
    prepared <- NULL
    if (!any(bad[names(prep_paths)])) {
      prepared_dir <- file.path(scratch, ftq_sanitize_label(subject))
      prep <- tryCatch(prepare_subject(subject, prepared_dir), error = identity)
      if (inherits(prep, "error")) {
        errors <- rbind(errors, failure("subject_input_preparation", conditionMessage(prep)))
      } else {
        prepared <- prep
        # Check the equivalent planned files in the actual destination.
        generated <- list.files(prepared_dir, recursive = TRUE)
        existing <- file.path(output_dir, ftq_sanitize_label(subject), "prepared", generated)
        if (length(generated)) {
          prep_error <- tryCatch(ftq_check_output_destinations(existing, overwrite), error = output_error)
          errors <- rbind(errors, prep_error)
        }
      }
    }
    for (focal in focal_regions) {
      checks <- if (!is.null(prepared)) tryCatch({
        manifest <- make_manifest(subject, focal, prepared)
        ftq_preflight(manifest)$checks
      }, error = function(e) failure("run_preflight", conditionMessage(e)))
      checks <- rbind(errors, checks)
      run_error <- tryCatch(check_run_destination(subject, focal), error = output_error)
      checks <- rbind(destination_error, run_error, checks)
      checks$detail <- gsub(scratch, "<temporary preparation>", checks$detail, fixed = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(subject = subject, focal_region = focal,
                                              checks, stringsAsFactors = FALSE)
    }
  }
  checks <- do.call(rbind, rows)
  rownames(checks) <- NULL
  index <- do.call(rbind, lapply(rows, function(x) {
    data.frame(subject = x$subject[1], focal_region = x$focal_region[1],
               passed = !any(x$severity == "error" & !x$passed), stringsAsFactors = FALSE)
  }))
  structure(list(passed = all(index$passed), index = index, checks = checks),
            class = "ftq_study_input_check")
}

#' @export
print.ftq_study_input_check <- function(x, ...) {
  failed <- unique(x$checks[!x$checks$passed & x$checks$severity == "error",
                            c("subject", "focal_region", "detail"), drop = FALSE])
  cat("Input checks passed for ", sum(x$index$passed), " of ", nrow(x$index),
      " requested analyses.\n", sep = "")
  if (nrow(failed)) print(failed, row.names = FALSE)
  notes <- sum(!x$checks$passed & x$checks$severity == "warning")
  if (notes) cat(notes, " non-blocking warning(s). See the checks table for details.\n", sep = "")
  cat("Review anatomical alignment independently before interpreting results.\n")
  invisible(x)
}
