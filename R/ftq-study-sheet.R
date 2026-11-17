# Convert one row per subject into the existing study input representation.
# File existence is checked by the workflow, so a check-only run can report
# missing files across subjects rather than dropping a row during import.
ftq_read_study_sheet <- function(study_sheet) {
  if (!is.character(study_sheet) || length(study_sheet) != 1L ||
      is.na(study_sheet) || !file.exists(study_sheet) || dir.exists(study_sheet)) {
    stop("study_sheet must point to an existing CSV file.", call. = FALSE)
  }
  sheet_dir <- dirname(normalizePath(study_sheet, mustWork = TRUE))
  lines <- readLines(study_sheet, encoding = "UTF-8", warn = FALSE)
  if (!length(lines) || !all(validUTF8(lines))) {
    stop("Save the study sheet as a UTF-8 CSV with a header and one row per subject.",
         call. = FALSE)
  }
  lines[1L] <- sub("^\ufeff", "", lines[1L])
  sheet <- tryCatch(utils::read.csv(
    text = lines, colClasses = "character", check.names = FALSE, row.names = NULL,
    fill = FALSE, strip.white = TRUE, na.strings = character()
  ), error = function(e) {
    stop("Could not read the study sheet ", study_sheet, ": ", conditionMessage(e),
         call. = FALSE)
  })
  headers <- names(sheet)
  if (any(!nzchar(headers)) || anyDuplicated(headers)) {
    stop("Each study sheet column needs a non-empty, unique header.", call. = FALSE)
  }
  missing <- setdiff(c("subject", "tractogram"), headers)
  if (length(missing)) {
    stop("Study sheet is missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  mask_columns <- headers[startsWith(headers, "mask:")]
  unknown <- setdiff(headers, c("subject", "tractogram", "reference_image", mask_columns))
  if (length(unknown)) {
    stop("Unrecognized study sheet columns: ", paste(unknown, collapse = ", "),
         ". Use subject, tractogram, optional reference_image, and mask:<region> columns.",
         call. = FALSE)
  }
  regions <- substring(mask_columns, 6L)
  if (length(regions) < 2L || any(!nzchar(trimws(regions)))) {
    stop("Supply at least two named region columns, such as mask:cerebellum and mask:cerebrum.",
         call. = FALSE)
  }
  ftq_regional_study_check_labels(regions, "Study sheet mask", reserved = "prepared")
  if (!nrow(sheet)) stop("The study sheet has no subject rows.", call. = FALSE)

  required <- sheet[c("subject", "tractogram", mask_columns)]
  empty <- which(do.call(cbind, lapply(required, function(x) !nzchar(trimws(x)))),
                 arr.ind = TRUE)
  if (length(empty)) {
    # Rows include the header so they can be located directly in a spreadsheet.
    cells <- paste0("row ", empty[, 1L] + 1L, " (", names(required)[empty[, 2L]], ")")
    stop("Fill in the required study sheet cells: ", paste(cells, collapse = ", "), ".",
         call. = FALSE)
  }
  subjects <- sheet$subject
  if (anyDuplicated(subjects)) {
    stop("Use one row per subject. Repeated subject IDs: ",
         paste(unique(subjects[duplicated(subjects)]), collapse = ", "), ".", call. = FALSE)
  }
  ftq_regional_study_check_labels(subjects, "Study sheet subject")
  path <- function(value) {
    # Accept relative paths copied from either Windows or macOS/Linux.
    ftq_manifest_path(gsub("\\", "/", value, fixed = TRUE), sheet_dir)
  }
  masks <- lapply(seq_len(nrow(sheet)), function(i) {
    stats::setNames(vapply(sheet[mask_columns], function(column) path(column[[i]]),
                           character(1L)), regions)
  })
  references <- NULL
  if ("reference_image" %in% headers && any(nzchar(trimws(sheet$reference_image)))) {
    references <- stats::setNames(lapply(sheet$reference_image, function(value) {
      if (nzchar(trimws(value))) path(value) else NULL
    }), subjects)
  }
  list(
    tractograms = stats::setNames(vapply(sheet$tractogram, path, character(1L)), subjects),
    region_masks = stats::setNames(masks, subjects),
    reference_images = references
  )
}
