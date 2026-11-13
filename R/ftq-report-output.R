# Reader-facing exports are separate from the canonical scientific result schema.
ftq_report_slices <- function(slices, nz = NULL) {
  if (is.null(slices)) return(NULL)
  if (!is.numeric(slices) || !length(slices) || any(!is.finite(slices)) ||
      any(slices < 1 | slices > .Machine$integer.max | slices != floor(slices)) || anyDuplicated(slices)) {
    stop("report_slices must contain distinct positive one-based axial indices.")
  }
  if (!is.null(nz) && any(slices > nz)) {
    stop("report_slices exceeds the image's axial dimension (", nz, ").")
  }
  as.integer(slices)
}

ftq_report_voxel_export <- function(result) {
  coords <- as.data.frame(ftq_coord_matrix(result$fNode_Coords))
  coordinate_names <- character()
  for (axis in c("x", "y", "z")) {
    name <- axis
    if (name %in% names(result)) name <- paste0("voxel_", axis)
    candidate <- name
    suffix <- 1L
    while (candidate %in% c(names(result), coordinate_names)) {
      suffix <- suffix + 1L
      candidate <- paste0(name, "_", suffix)
    }
    coordinate_names <- c(coordinate_names, candidate)
  }
  names(coords) <- coordinate_names
  reserved <- c("fNode_Coords", "tracks", "local", "leave_no_target", "none", "cNode_Coords")
  roles <- ifelse(names(result) %in% c("fNode_Coords", "cNode_Coords"), "identifier",
                  ifelse(names(result) %in% reserved, "count", "target"))
  columns <- data.frame(
    column = c(coordinate_names, names(result)),
    role = c(rep("coordinate", 3L), roles),
    source_column = c(rep("fNode_Coords", 3L), names(result)),
    units = c(rep("one-based focal-image array index", 3L),
              ifelse(roles == "identifier", "coordinate key", "assigned streamline count")),
    description = c(paste("Focal-image", c("x", "y", "z"), "array index; not scanner millimetres."),
                    paste("Unchanged canonical field:", names(result))),
    stringsAsFactors = FALSE
  )
  columns$units[columns$column == "none" & columns$role == "count"] <- "binary zero-assignment indicator"
  list(data = cbind(coords, result), columns = columns)
}

ftq_report_local_fraction <- function(result) {
  positive <- is.finite(result$tracks) & result$tracks > 0
  fraction <- rep(NA_real_, nrow(result))
  fraction[positive] <- result$local[positive] / result$tracks[positive]
  if (any(!is.finite(fraction[positive]) | fraction[positive] < 0 | fraction[positive] > 1)) {
    stop("Local fractions require finite local counts between zero and tracks.")
  }
  data.frame(fNode_Coords = result$fNode_Coords, local_fraction = fraction,
             has_assigned_streamlines = positive, stringsAsFactors = FALSE)
}

ftq_report_write_voxel_exports <- function(result, output_dir, prefix,
                                           template_image = NULL, overwrite = FALSE) {
  table_dir <- file.path(output_dir, "tables")
  report_paths <- ftq_report_output_paths(output_dir, prefix)
  tables <- report_paths$tables[c("all_voxels", "all_voxel_columns", "local_fraction")]
  maps <- list()
  if (!is.null(template_image) && file.exists(template_image)) {
    maps <- report_paths$maps
  }
  ftq_check_overwrite(c(unlist(tables), unlist(maps)), overwrite)
  view <- ftq_report_voxel_export(result)
  fraction <- ftq_report_local_fraction(result)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(view$data, tables$all_voxels, row.names = FALSE)
  utils::write.csv(view$columns, tables$all_voxel_columns, row.names = FALSE)
  utils::write.csv(fraction, tables$local_fraction, row.names = FALSE, na = "NA")
  if (length(maps)) {
    template <- RNifti::readNifti(template_image, internal = FALSE)
    dims <- dim(template)[1:3]
    coords <- ftq_coord_matrix(result$fNode_Coords)
    if (any(coords < 1) || any(sweep(coords, 2L, dims, ">"))) {
      stop("Focal coordinates are outside the local-fraction template grid.")
    }
    dir.create(dirname(maps$local_fraction), recursive = TRUE, showWarnings = FALSE)
    values <- array(0, dim = dims)
    valid <- array(0L, dim = dims)
    values[coords] <- ifelse(fraction$has_assigned_streamlines, fraction$local_fraction, 0)
    valid[coords] <- as.integer(fraction$has_assigned_streamlines)
    RNifti::writeNifti(values, maps$local_fraction, template = template, datatype = "float")
    RNifti::writeNifti(valid, maps$local_fraction_valid, template = template, datatype = "uint8")
  }
  list(tables = tables, maps = maps)
}

# Explicit roles, rather than a new reserved-name blacklist, allow genuine
# targets named x/y/z to survive a reader-facing CSV round trip.
ftq_report_read_voxel_view <- function(path, table) {
  if (!grepl("_all-voxels[.]csv$", path)) return(table)
  schema_path <- sub("_all-voxels[.]csv$", "_all-voxels-columns.csv", path)
  if (!file.exists(schema_path)) {
    stop("The all-voxel view needs its companion all-voxels-columns.csv. Use the canonical per-voxel-hits.csv for analysis.")
  }
  schema <- utils::read.csv(schema_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("column", "role", "source_column") %in% names(schema)) ||
      anyDuplicated(schema$column) || !identical(as.character(schema$column), names(table)) ||
      any(!schema$role %in% c("coordinate", "identifier", "count", "target"))) {
    stop("The all-voxel column-role record does not match its CSV.")
  }
  keep <- schema$role != "coordinate"
  if (anyDuplicated(schema$source_column[keep])) stop("Ambiguous source fields in the all-voxel column-role record.")
  result <- table[schema$column[keep]]
  names(result) <- schema$source_column[keep]
  if (!"fNode_Coords" %in% names(result) || sum(!keep) != 3L) stop("Incomplete all-voxel coordinate record.")
  coords <- ftq_coord_matrix(result$fNode_Coords)
  if (!isTRUE(all.equal(unname(as.matrix(table[schema$column[!keep]])),
                       unname(coords), check.attributes = FALSE))) {
    stop("The all-voxel coordinates do not reconstruct fNode_Coords.")
  }
  result
}

ftq_report_html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

ftq_report_downloads <- function(paths, root, labels = names(paths)) {
  if (!length(paths)) return("")
  items <- vapply(seq_along(paths), function(i) {
    p <- paths[[i]]
    if (!is.character(p) || length(p) != 1L || is.na(p) || !file.exists(p)) return("")
    paste0("<li><a href='", ftq_report_html_escape(ftq_report_href(p, root)), "'>",
           ftq_report_html_escape(labels[[i]]), "</a></li>")
  }, character(1))
  if (!any(nzchar(items))) "" else paste0("<ul>", paste(items, collapse = ""), "</ul>")
}

ftq_report_metadata_html <- function(metadata_path, root) {
  if (is.null(metadata_path) || !file.exists(metadata_path)) return("")
  parsed <- tryCatch(jsonlite::read_json(metadata_path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(parsed)) return("<p class='warn'>Saved metadata could not be read. The original JSON is available in downloads.</p>")
  rows <- list()
  walk <- function(value, key) {
    if (is.list(value) && length(value)) {
      keys <- names(value)
      if (is.null(keys)) keys <- as.character(seq_along(value))
      for (i in seq_along(value)) walk(value[[i]], paste(key, keys[[i]], sep = " / "))
    } else {
      rows[[length(rows) + 1L]] <<- paste0("<tr><th>", ftq_report_html_escape(gsub("_", " ", key, fixed = TRUE)),
                                         "</th><td>", ftq_report_html_escape(if (length(value)) paste(value, collapse = ", ") else "Not supplied"), "</td></tr>")
    }
  }
  for (key in intersect(c("generated_at", "package", "manifest", "inputs", "summary", "endpoint_analysis", "local_search_audit"), names(parsed))) {
    walk(parsed[[key]], key)
  }
  context <- parsed$files$context
  if (is.list(context)) {
    for (key in intersect(c("source", "summary", "metadata"), names(context))) {
      walk(context[[key]], paste("tractogram context", key, sep = " / "))
    }
  }
  paste0("<div class='table-scroll' role='region' aria-label='Saved run provenance' tabindex='0'><table>",
         paste(rows, collapse = ""), "</table></div>",
         ftq_report_downloads(list(metadata_path), root, "Complete run metadata (JSON)"))
}
