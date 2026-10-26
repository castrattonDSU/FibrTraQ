# Resolve report labels from saved evidence before the ordinary manifest
# normalizer supplies historical defaults. This reads metadata only; it never
# imports result tables, maps, endpoint audits, or other analysis files.
ftq_report_resolve_manifest <- function(result_csv, manifest = NULL) {
  scalar <- function(x) {
    if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(as.character(x))) return(NULL)
    if (is.character(x) && toupper(trimws(x)) %in% c("", "NA", "NULL")) return(NULL)
    x
  }
  explicit <- manifest
  if (!is.null(explicit) && !is.list(explicit)) {
    if (!is.character(explicit) || length(explicit) != 1L) stop("manifest must be a path or list.")
    path <- path.expand(explicit)
    if (!file.exists(path)) stop("Manifest file does not exist: ", explicit)
    explicit <- if (grepl("[.]json$", path, ignore.case = TRUE)) {
      jsonlite::read_json(path, simplifyVector = FALSE)
    } else if (grepl("[.]ya?ml$", path, ignore.case = TRUE)) {
      yaml::read_yaml(path)
    } else stop("Manifest must be a .json, .yml, or .yaml file.")
    # Preserve the existing path-manifest interface and path normalization.
    manifest <- ftq_read_manifest(manifest)
  }
  if (!is.null(explicit) && !is.list(explicit)) stop("Manifest content must be a list.")

  metadata <- NULL
  metadata_path <- NULL
  if (is.character(result_csv) && length(result_csv) == 1L &&
      grepl("_fibrtraq_per-voxel-hits[.]csv$", result_csv)) {
    candidate <- sub("_fibrtraq_per-voxel-hits[.]csv$", "_fibrtraq_run-metadata.json", result_csv)
    if (file.exists(candidate)) {
      metadata <- tryCatch(jsonlite::read_json(candidate, simplifyVector = FALSE),
        error = function(e) stop("Cannot read companion run metadata: ", conditionMessage(e)))
      if (!is.list(metadata)) stop("Companion run metadata must contain a JSON object.")
      recorded_file <- scalar(metadata$files$result_csv)
      if (!is.null(recorded_file) && !identical(
        basename(ftq_report_normalize_path(recorded_file)), basename(result_csv)
      )) stop("Companion run metadata names a different result CSV.")
      metadata_path <- candidate
    }
  }
  saved <- metadata$manifest
  if (!is.null(saved) && !is.list(saved)) stop("Companion run manifest must be a JSON object.")
  result_value <- function(field) {
    if (!is.data.frame(result_csv)) return(NULL)
    scalar(attr(result_csv, paste0("fibrtraq_", field), exact = TRUE))
  }
  manifest_value <- function(x, field) ftq_value(scalar(x$analysis[[field]]), scalar(x[[field]]))
  recorded_value <- function(field, normalize = identity, additional = NULL) {
    values <- Filter(Negate(is.null), list(result_value(field), scalar(metadata$summary[[field]]),
      scalar(saved$analysis[[field]]), scalar(saved[[field]]), scalar(additional)))
    values <- lapply(values, normalize)
    if (length(unique(values)) > 1L) stop("Conflicting ", field, " values in recorded result provenance.")
    if (length(values)) values[[1L]] else NULL
  }
  recorded_mode <- recorded_value("target_intersection_mode", ftq_normalize_target_intersection_mode)
  supplied_mode <- manifest_value(explicit, "target_intersection_mode")
  if (!is.null(supplied_mode)) supplied_mode <- ftq_normalize_target_intersection_mode(supplied_mode)
  if (!is.null(supplied_mode) && !is.null(recorded_mode) && !identical(supplied_mode, recorded_mode)) {
    stop("Manifest target_intersection_mode conflicts with recorded result provenance (", recorded_mode, ").")
  }
  mode <- ftq_value(supplied_mode, recorded_mode)
  notes <- character()
  if (is.null(mode)) {
    mode <- "full_streamline"
    notes <- c(notes, paste0("Target-assignment method was not recorded. full_streamline is the historical ",
      "report default, not a verified method for these data. Supply the original ",
      "target_intersection_mode or retain the matching run-metadata JSON. ",
      "Endpoint QC settings do not establish the target-assignment method."))
  }
  supplied_radius <- manifest_value(explicit, "local_search_radius_mm")
  recorded_radius <- if (identical(mode, "local_search")) {
    recorded_value("local_search_radius_mm", ftq_validate_local_search_radius,
      metadata$local_search_audit$search_radius_mm)
  } else NULL
  if (!is.null(supplied_radius)) supplied_radius <- ftq_validate_local_search_radius(supplied_radius)
  if (identical(mode, "local_search") && !is.null(supplied_radius) && !is.null(recorded_radius) &&
      !isTRUE(all.equal(supplied_radius, recorded_radius))) {
    stop("Manifest local_search_radius_mm conflicts with recorded result provenance (", recorded_radius, " mm).")
  }
  radius <- ftq_value(supplied_radius, recorded_radius)
  if (is.null(radius)) {
    radius <- 2
    if (identical(mode, "local_search")) notes <- c(notes,
      "Local-search radius was not recorded. The historical 2 mm report default is unverified; supply the original local_search_radius_mm.")
  }
  method_choices <- list(
    focal_assignment_mode = c("all_focal_voxels", "first_intersection"),
    streamline_voxelization_mode = c("sample_points", "segment_interpolation")
  )
  other_methods <- list()
  for (field in names(method_choices)) {
    normalize <- function(value) match.arg(as.character(value), method_choices[[field]])
    recorded <- recorded_value(field, normalize)
    supplied <- manifest_value(explicit, field)
    if (!is.null(supplied)) supplied <- normalize(supplied)
    if (!is.null(supplied) && !is.null(recorded) && !identical(supplied, recorded)) {
      stop("Manifest ", field, " conflicts with recorded result provenance (", recorded, ").")
    }
    other_methods[[field]] <- ftq_value(supplied, recorded)
  }
  # Supplied identity and presentation settings take precedence. Fill missing
  # fields from saved provenance without guessing paths in a moved bundle.
  resolved <- if (is.null(manifest)) ftq_value(saved, list()) else manifest
  for (field in setdiff(names(saved), "analysis")) {
    if (is.null(resolved[[field]])) resolved[[field]] <- saved[[field]]
  }
  if (!is.null(manifest)) {
    for (field in names(saved$analysis)) {
      if (is.null(manifest_value(explicit, field))) {
        resolved$analysis[[field]] <- saved$analysis[[field]]
      }
    }
  }
  for (field in names(method_choices)) {
    value <- other_methods[[field]]
    # Drop missing sentinels from saved manifests as well as summary records;
    # ordinary report normalization can then apply the legacy default.
    resolved[[field]] <- value
    if (!is.null(resolved$analysis)) resolved$analysis[[field]] <- value
  }
  resolved$target_intersection_mode <- mode
  resolved$local_search_radius_mm <- radius
  # Nested analysis fields otherwise take priority in ftq_report's normalizer.
  if (!is.null(resolved$analysis)) {
    resolved$analysis$target_intersection_mode <- mode
    resolved$analysis$local_search_radius_mm <- radius
  }
  resolved$provenance_warning <- if (length(notes)) paste(notes, collapse = " ") else NULL
  resolved$provenance_metadata_json <- metadata_path
  if (length(notes) && !is.data.frame(result_csv)) warning(resolved$provenance_warning, call. = FALSE)
  resolved
}
