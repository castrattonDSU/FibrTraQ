#!/usr/bin/env Rscript
# Reproduce the native-TCK public-data validation without changing its inputs.
# Rscript fibrtraq-tck-public-validation.R /path/to/prepared-inputs /path/to/new-results
# Install the tested FibrTraQ revision first. Acquisition and mask-preparation
# provenance accompany the prepared inputs. This runner does not download data.

args <- commandArgs(trailingOnly = TRUE)
oracle_only <- "--oracle-only" %in% args
args <- args[args != "--oracle-only"]
if (length(args) != 2L || any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript fibrtraq-tck-public-validation.R <input_dir> <new_output_dir> [--oracle-only]\n")
  quit(status = if (any(args %in% c("--help", "-h"))) 0L else 1L)
}
suppressPackageStartupMessages(library(FibrTraQ))
if (!exists("ftq_read_tractogram", envir = asNamespace("FibrTraQ"), inherits = FALSE)) {
  stop("Install the FibrTraQ revision with native TCK support first.", call. = FALSE)
}
input_dir <- normalizePath(path.expand(args[[1L]]), mustWork = TRUE)
output_dir <- normalizePath(path.expand(args[[2L]]), mustWork = FALSE)
if (!dir.exists(input_dir)) stop("input_dir must be a directory.", call. = FALSE)
if (file.exists(output_dir)) {
  stop("output_dir already exists. Choose a NEW directory to preserve prior outputs.", call. = FALSE)
}
inputs <- c(
  tck = "tracks_10k_6fed19354195.tck", trk = "native-equivalent.trk",
  reference = "reference.nii.gz", focal = "focal-left-interface.nii.gz",
  positive_x_cortical_GM = "target-right-gm.nii.gz",
  negative_x_cortical_GM = "target-left-gm.nii.gz",
  reader_oracle = "reader-oracle.json"
)
if (oracle_only) inputs <- inputs[c("tck", "reader_oracle")]
inputs <- stats::setNames(file.path(input_dir, inputs), names(inputs))
missing <- inputs[!file.exists(inputs) | file.info(inputs)$isdir]
if (length(missing)) stop("Missing prepared input files: ", paste(missing, collapse = ", "))
inputs <- vapply(inputs, normalizePath, character(1), mustWork = TRUE)
if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create output_dir: ", output_dir, call. = FALSE)
}
# normalizePath() can leave a nonexistent relative directory unchanged. Resolve
# it after creation so reused manifests cannot nest outputs under dataset_root.
output_dir <- normalizePath(output_dir, mustWork = TRUE)

# Counts, RSC values, CSV tables, and density/map values must match exactly.
# Only continuous rendered world coordinates may differ by 0.0001 mm after
# equivalent TRK Float32 serialization; no tolerance applies to counts.
viewer_tolerance_mm <- 1e-4
image_affine_tolerance <- 1e-6
format_provenance_columns <- c("selected_track_space", "track_voxel_size_mm",
                               "mapping_selection_rule", "backend", "reverse_streamlines")
parity <- data.frame(comparison = character(), artifact = character(),
                     passed = logical(), tolerance = numeric(),
                     max_absolute_difference = numeric(), detail = character())
state <- list(
  status = "running", started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  validation_mode = if (oracle_only) "independent_reader_oracle_only" else "complete_public_validation",
  dataset = "Public Neurodesk sub-CON02 10,000-streamline TCK example",
  scope = paste("One public human dataset and one biological subject; repeated",
                "computational runs are not additional biological samples.",
                "The masks are a priori tissue/RAS-x geometric demonstrations,",
                "not validated hemispheric labels or atlas parcels.",
                "This validates software interoperability, not connectivity ground truth."),
  input_directory = input_dir,
  inputs = data.frame(name = names(inputs), filename = basename(inputs),
                      bytes = as.numeric(file.info(inputs)$size),
                      md5 = unname(tools::md5sum(inputs))),
  package = as.list(utils::packageDescription("FibrTraQ")[c(
    "Package", "Version", "RemoteRepo", "RemoteRef", "RemoteSha", "Built"
  )]),
  tolerances = list(counts_and_tables = 0, nifti_voxel_values = 0,
                    nifti_affine_matrix = image_affine_tolerance,
                    viewer_world_coordinates_mm = viewer_tolerance_mm),
  excluded_roi_summary_columns = format_provenance_columns,
  excluded_columns_reason = "Coordinate interpretation, backend, and requested reversal differ by design; scientific values remain exact.",
  runs = list(), warnings = character(), checks = parity
)
summarize_validation <- function() {
  native_comparisons <- c(
    "independent_reader_oracle", "independent_reader_oracle_stream",
    "independent_reader_oracle_memory", "tck_stream_cpp", "public_input",
    "tck-memory-r", "tck-reversed-stream-cpp"
  )
  format_comparisons <- c("trk_stream_cpp", "native_tck_vs_equivalent_trk", "full_report_exports")
  summarize <- function(rows, required_comparisons) {
    missing <- setdiff(required_comparisons, rows$comparison)
    failed <- rows[!rows$passed, c("comparison", "artifact"), drop = FALSE]
    list(
      status = if (!nrow(rows)) "not_run" else if (nrow(failed)) "failed" else
        if (length(missing)) "incomplete" else "passed",
      checks_run = nrow(rows), checks_passed = sum(rows$passed),
      checks_failed = nrow(failed), missing_comparisons = missing,
      failed_artifacts = failed
    )
  }
  native_rows <- parity[parity$comparison %in% native_comparisons |
    (parity$comparison == "execution_provenance" &
       parity$artifact == "native_backends_and_reversal"), , drop = FALSE]
  format_rows <- parity[parity$comparison %in% format_comparisons |
    (parity$comparison == "execution_provenance" &
       parity$artifact == "equivalent_trk_coordinate_space"), , drop = FALSE]
  list(
    native_tck_execution_and_invariance = summarize(
      native_rows, c(native_comparisons, "execution_provenance")
    ),
    tck_vs_trk_representation_equivalence = summarize(
      format_rows, c(format_comparisons, "execution_provenance")
    )
  )
}
write_state <- function() {
  state$checks <- parity
  state$validation_summary <- summarize_validation()
  jsonlite::write_json(state, file.path(output_dir, "validation-results.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  utils::write.csv(parity, file.path(output_dir, "parity-checks.csv"), row.names = FALSE)
}
record <- function(comparison, artifact, passed, tolerance = 0,
                   max_difference = NA_real_, detail = "") {
  parity <<- rbind(parity, data.frame(
    comparison = comparison, artifact = artifact, passed = isTRUE(passed),
    tolerance = tolerance, max_absolute_difference = max_difference,
    detail = detail, stringsAsFactors = FALSE
  ))
  write_state()
}
read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
check_table <- function(a, b, comparison, artifact, exclude = character()) {
  a <- read_csv(a)
  b <- read_csv(b)
  a[intersect(exclude, names(a))] <- NULL
  b[intersect(exclude, names(b))] <- NULL
  difference <- all.equal(a, b, tolerance = 0, check.attributes = TRUE)
  numeric_columns <- intersect(names(a)[vapply(a, is.numeric, logical(1))],
                               names(b)[vapply(b, is.numeric, logical(1))])
  maximum <- NA_real_
  if (identical(dim(a), dim(b)) && length(numeric_columns)) {
    delta <- abs(unlist(a[numeric_columns]) - unlist(b[numeric_columns]))
    if (any(is.finite(delta))) maximum <- max(delta[is.finite(delta)])
  }
  record(comparison, artifact, isTRUE(difference), max_difference = maximum,
         detail = if (isTRUE(difference)) "Exact values, row order, labels and denominators match." else
           paste(difference, collapse = "; "))
}
check_image <- function(a, b, comparison, artifact) {
  a <- RNifti::readNifti(a, internal = FALSE)
  b <- RNifti::readNifti(b, internal = FALSE)
  dims_equal <- identical(dim(a), dim(b))
  values_equal <- dims_equal && identical(as.numeric(a), as.numeric(b))
  maximum <- if (dims_equal) max(abs(as.numeric(a) - as.numeric(b))) else NA_real_
  record(comparison, paste0(artifact, ":voxel_values"), values_equal,
         max_difference = maximum)
  affine_difference <- max(abs(as.numeric(RNifti::xform(a)) - as.numeric(RNifti::xform(b))))
  record(comparison, paste0(artifact, ":grid"),
         dims_equal && is.finite(affine_difference) &&
           affine_difference <= image_affine_tolerance &&
           identical(RNifti::pixunits(a), RNifti::pixunits(b)),
         tolerance = image_affine_tolerance, max_difference = affine_difference)
}
read_scene <- function(path) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  pattern <- '<script id="scene-data" type="application/json">(.*?)</script>'
  matched <- regmatches(html, regexec(pattern, html, perl = TRUE))[[1L]]
  if (length(matched) != 2L) stop("Missing 3D viewer scene data: ", path)
  jsonlite::fromJSON(matched[[2L]], simplifyVector = FALSE)
}
# Strings, labels, logical values, list order and list lengths must match.
# Numeric scene values use an absolute bound; integer count differences of one
# or more also fail because the bound is less than one.
check_scene <- function(a, b, comparison) {
  maximum <- 0
  mismatches <- character()
  compare <- function(x, y, location = "scene") {
    if (is.numeric(x) && is.numeric(y) && length(x) == length(y)) {
      if (!identical(is.na(x), is.na(y))) {
        mismatches <<- c(mismatches, location)
      } else {
        differences <- abs(x - y)
        if (any(!is.finite(differences[!is.na(differences)]))) {
          mismatches <<- c(mismatches, location)
        } else if (any(!is.na(differences))) {
          maximum <<- max(maximum, differences, na.rm = TRUE)
          if (any(differences > viewer_tolerance_mm, na.rm = TRUE)) {
            mismatches <<- c(mismatches, location)
          }
        }
      }
    } else if (is.list(x) && is.list(y) &&
               identical(names(x), names(y)) && length(x) == length(y)) {
      for (i in seq_along(x)) compare(x[[i]], y[[i]], paste0(location, "/", i))
    } else if (!identical(x, y)) {
      mismatches <<- c(mismatches, location)
    }
  }
  compare(read_scene(a), read_scene(b))
  record(comparison, "3d_viewer_geometry_and_annotations", !length(mismatches),
         tolerance = viewer_tolerance_mm, max_difference = maximum,
         detail = if (length(mismatches)) paste(head(mismatches, 30L), collapse = "; ") else
           "Continuous geometry within stated absolute tolerance; annotations and structure match.")
}
time_run <- function(label, expression) {
  message("Running ", label)
  start <- proc.time()[["elapsed"]]
  value <- force(expression)
  state$runs[[label]] <<- list(elapsed_seconds = proc.time()[["elapsed"]] - start,
                              status = "completed")
  write_state()
  value
}

validate_reader_oracle <- function() {
  oracle <- read_json(inputs[["reader_oracle"]])
  coordinate_encoding <- "Float64LE; C-order xyz triples, concatenated in stored streamline/point order; no delimiters"
  point_count_encoding <- "Int32LE; one point count per stored streamline"
  digest_valid <- function(x, digits) {
    is.character(x) && length(x) == 1L && !is.na(x) &&
      grepl(paste0("^[[:xdigit:]]{", digits, "}$"), x)
  }
  declared <- identical(as.integer(oracle$schema_version), 1L) &&
    identical(as.numeric(oracle$streamlines), 10000) &&
    identical(as.numeric(oracle$points), 272137) &&
    identical(oracle$coordinate_encoding, coordinate_encoding) &&
    identical(oracle$point_count_encoding, point_count_encoding) &&
    identical(as.numeric(oracle$coordinates$bytes), 272137 * 3 * 8) &&
    identical(as.numeric(oracle$points_per_streamline$bytes), 10000 * 4) &&
    digest_valid(oracle$coordinates$md5, 32L) &&
    digest_valid(oracle$coordinates$sha256, 64L) &&
    digest_valid(oracle$points_per_streamline$md5, 32L) &&
    digest_valid(oracle$points_per_streamline$sha256, 64L)
  record("independent_reader_oracle", "declaration", declared,
         detail = "NiBabel oracle: 10,000 streamlines, 272,137 points; canonical binary encodings and digest lengths validated.")
  if (!declared) stop("Invalid or unexpected independent reader oracle; no analyses were run.")
  state$reader_oracle <<- list(reference = oracle, actual = list())
  read_tracks <- getFromNamespace("ftq_read_tractogram", "FibrTraQ")
  iterate_tracks <- getFromNamespace("ftq_track_stream_iterate", "FibrTraQ")
  canonical_digest <- function(reader) {
    coordinates_file <- tempfile("tck-oracle-coordinates-", fileext = ".bin")
    point_counts_file <- tempfile("tck-oracle-point-counts-", fileext = ".bin")
    coordinates_con <- file(coordinates_file, "wb")
    point_counts_con <- file(point_counts_file, "wb")
    on.exit({
      if (!is.null(coordinates_con)) close(coordinates_con)
      if (!is.null(point_counts_con)) close(point_counts_con)
      unlink(c(coordinates_file, point_counts_file))
    }, add = TRUE)
    streamline_count <- 0
    point_count <- 0
    append_tracks <- function(chunk) {
      for (track in chunk) {
        coords <- as.matrix(track$coords)
        if (ncol(coords) != 3L || any(!is.finite(coords))) {
          stop("Reader oracle encountered invalid point coordinates.")
        }
        writeBin(as.numeric(t(coords)), coordinates_con, size = 8L, endian = "little")
        writeBin(as.integer(nrow(coords)), point_counts_con, size = 4L, endian = "little")
        streamline_count <<- streamline_count + 1
        point_count <<- point_count + nrow(coords)
      }
      invisible(NULL)
    }
    tracks <- read_tracks(inputs[["tck"]], reader = reader, chunk_size = 997L)
    if (reader == "stream") {
      iterate_tracks(tracks, function(chunk, start_index, end_index) append_tracks(chunk))
    } else {
      append_tracks(tracks$tracks)
    }
    close(coordinates_con)
    coordinates_con <- NULL
    close(point_counts_con)
    point_counts_con <- NULL
    list(
      streamlines = streamline_count, points = point_count,
      coordinates = list(bytes = as.numeric(file.info(coordinates_file)$size),
                         md5 = unname(tools::md5sum(coordinates_file))),
      points_per_streamline = list(bytes = as.numeric(file.info(point_counts_file)$size),
                                   md5 = unname(tools::md5sum(point_counts_file)))
    )
  }
  for (reader in c("stream", "memory")) {
    actual <- time_run(paste0("independent_reader_oracle_", reader), canonical_digest(reader))
    state$reader_oracle$actual[[reader]] <<- actual
    label <- paste0("independent_reader_oracle_", reader)
    record(label, "streamline_count", identical(actual$streamlines, 10000),
           max_difference = abs(actual$streamlines - 10000))
    record(label, "point_count", identical(actual$points, 272137),
           max_difference = abs(actual$points - 272137))
    for (component in c("coordinates", "points_per_streamline")) {
      record(label, paste0(component, "_byte_exact"),
             identical(actual[[component]]$bytes, as.numeric(oracle[[component]]$bytes)) &&
               identical(actual[[component]]$md5, tolower(oracle[[component]]$md5)),
             detail = paste0("R MD5=", actual[[component]]$md5,
                             "; independent NiBabel MD5=", oracle[[component]]$md5,
                             "; bytes=", actual[[component]]$bytes))
    }
  }
  passed <- all(parity$passed[grepl("^independent_reader_oracle", parity$comparison)])
  jsonlite::write_json(list(
    passed = passed, comparison = "R native TCK readers versus independently serialized NiBabel coordinates",
    tolerance = 0, coordinate_encoding = coordinate_encoding,
    point_count_encoding = point_count_encoding, oracle = state$reader_oracle,
    scope = "All stored coordinates and streamline lengths, in original order; this is a parser check, not an analysis or report check."
  ), file.path(output_dir, "reader-oracle-results.json"),
  auto_unbox = TRUE, pretty = TRUE, null = "null")
  if (!passed) stop("Independent reader oracle mismatch; no analyses were run.")
  invisible(TRUE)
}

main <- function() {
  validate_reader_oracle()
  if (oracle_only) return(invisible(NULL))
  set.seed(20261112)
  targets <- inputs[c("positive_x_cortical_GM", "negative_x_cortical_GM")]
  make_full <- function(format) ftq_run_files(
    tractogram = inputs[[format]], focal_mask = inputs[["focal"]],
    target_masks = targets, reference_image = inputs[["reference"]],
    output_dir = file.path(output_dir, paste0(format, "-stream-cpp")),
    subject = "sub-CON02", focal_label = "negative_x_interface", description = "public-tck-validation",
    track_space = if (format == "tck") "auto" else "trackvis",
    axis_flip = if (format == "tck") NULL else rep(FALSE, 3L),
    discretization = "round", assignment_mode = "full_streamline",
    focal_assignment_mode = "all_focal_voxels",
    streamline_voxelization_mode = "segment_interpolation",
    endpoint_analysis = "exact_and_dilated1", track_reader = "stream",
    track_chunk_size = 1000L, backend = "cpp", sample_per_target = 25L,
    overwrite = FALSE, verbose = TRUE, progress = TRUE, open = FALSE
  )
  native <- time_run("tck_stream_cpp", make_full("tck"))
  equivalent <- time_run("trk_stream_cpp", make_full("trk"))
  for (bundle in list(native, equivalent)) {
    label <- if (identical(bundle$output_dir, native$output_dir)) "tck_stream_cpp" else "trk_stream_cpp"
    record(label, "complete_report_bundle", nrow(bundle$index) == 1L &&
             identical(bundle$index$run_status, "completed") &&
             isTRUE(bundle$index$validation_passed) &&
             all(file.exists(c(bundle$report_index, bundle$index$primary_report,
                               bundle$index$exploratory_report, bundle$index$interactive_3d_qc))))
  }
  metadata <- read_json(native$index$metadata_json)
  metadata_trk <- read_json(equivalent$index$metadata_json)
  baseline <- as.list(native$index[1L, ])
  compare_outputs <- function(actual, maps, label) {
    keys <- c("result_csv", "rsc_csv", "fate_csv", "roi_summary_csv",
              "endpoint_summary_csv", "endpoint_penetration_csv", "endpoint_context_csv")
    for (key in keys) {
      check_table(baseline[[key]], actual[[key]], label, key,
                  exclude = if (key == "roi_summary_csv") format_provenance_columns else character())
    }
    expected_maps <- metadata$files$map_paths
    record(label, "map_names", setequal(names(expected_maps), names(maps)))
    for (name in intersect(names(expected_maps), names(maps))) {
      check_image(expected_maps[[name]], maps[[name]], label, paste0("map:", name))
    }
  }
  compare_outputs(as.list(equivalent$index[1L, ]), metadata_trk$files$map_paths,
                  "native_tck_vs_equivalent_trk")
  check_scene(native$index$interactive_3d_qc, equivalent$index$interactive_3d_qc,
              "native_tck_vs_equivalent_trk")
  # Full reports use supplied anatomy. Generate density explicitly through the
  # public API so the alternative tractogram-derived context is validated too.
  contexts <- lapply(c("tck", "trk"), function(format) time_run(
    paste0(format, "_density_context"), ftq_tractogram_context(
      tractogram = inputs[[format]], template_image = inputs[["focal"]],
      ref_image = inputs[["reference"]],
      output_dir = file.path(output_dir, paste0(format, "-density-context")),
      prefix = "sub-CON02", track_space = if (format == "tck") "auto" else "trackvis",
      axis_flip = rep(FALSE, 3L), track_chunk_size = 1000L,
      density_discretization = "round", streamline_voxelization_mode = "segment_interpolation",
      overwrite = FALSE
    )
  ))
  check_image(contexts[[1L]]$density_image, contexts[[2L]]$density_image,
              "native_tck_vs_equivalent_trk", "tractogram_density")
  check_table(contexts[[1L]]$summary_csv, contexts[[2L]]$summary_csv,
              "native_tck_vs_equivalent_trk", "tractogram_density_summary")
  record("public_input", "processed_10000_streamlines",
         identical(as.numeric(contexts[[1L]]$summary$processed_streamlines), 10000))

  # Reuse the exact normalized native-run manifest for lighter checks. Every
  # scientific output remains enabled; only redundant reports are omitted.
  light_run <- function(label, reader, backend, reverse = FALSE) {
    manifest <- metadata$manifest
    manifest$execution <- list(track_reader = reader, track_chunk_size = 1000L)
    manifest$analysis$reverse_streamlines <- reverse
    manifest$outputs <- list(output_dir = file.path(output_dir, label),
                             write_maps = TRUE, write_report = FALSE, overwrite = FALSE)
    manifest$report$context_source <- "none"
    run <- time_run(label, ftq_run(manifest, backend = backend, verbose = TRUE, progress = TRUE))
    record(label, "computational_validation", isTRUE(run$validation$passed))
    compare_outputs(run$outputs, run$outputs$map_paths, label)
    run
  }
  memory <- light_run("tck-memory-r", "memory", "r")
  reversed <- light_run("tck-reversed-stream-cpp", "stream", "cpp", reverse = TRUE)
  native_roi <- read_csv(baseline$roi_summary_csv)
  memory_roi <- read_csv(memory$outputs$roi_summary_csv)
  reversed_roi <- read_csv(reversed$outputs$roi_summary_csv)
  record("execution_provenance", "native_backends_and_reversal",
         identical(native_roi$backend, "cpp") && identical(memory_roi$backend, "r") &&
           identical(reversed_roi$backend, "cpp") && isTRUE(reversed_roi$reverse_streamlines) &&
           identical(native_roi$selected_track_space, "world"))
  record("execution_provenance", "equivalent_trk_coordinate_space",
         identical(read_csv(equivalent$index$roi_summary_csv)$selected_track_space, "trackvis"))
  # Compare the additional full-report tables, including the target LUT.
  a_dir <- dirname(native$index$result_csv)
  b_dir <- dirname(equivalent$index$result_csv)
  a_csv <- list.files(a_dir, "[.]csv$", recursive = TRUE)
  b_csv <- list.files(b_dir, "[.]csv$", recursive = TRUE)
  record("full_report_exports", "matching_csv_file_set", setequal(a_csv, b_csv))
  already_compared <- basename(unlist(baseline[c("result_csv", "rsc_csv", "fate_csv",
    "roi_summary_csv", "endpoint_summary_csv", "endpoint_penetration_csv", "endpoint_context_csv")]))
  for (path in intersect(a_csv, b_csv)) {
    if (!basename(path) %in% already_compared) {
      check_table(file.path(a_dir, path), file.path(b_dir, path), "full_report_exports", path)
    }
  }
  science <- list(
    dataset = state$dataset, biological_subjects = 1L, independent_datasets = 1L,
    primary_computational_runs = 4L, density_context_runs = 2L,
    masks = list(
      negative_x_interface = "GM/WM interface >0.5 and RAS world x<0",
      positive_x_cortical_GM = "Cortical GM >0.5 and RAS world x>=0",
      negative_x_cortical_GM = "Cortical GM >0.5 and RAS world x<0, excluding the focal mask",
      interpretation = "A priori tissue/geometric demonstration masks, not validated hemispheres or atlas parcels."),
    analysis = list(assignment_mode = "full_streamline", focal_assignment_mode = "all_focal_voxels",
                     streamline_voxelization_mode = "segment_interpolation", discretization = "round",
                     endpoint_analysis = "exact_and_dilated1"),
    summary = metadata$summary,
    roi_summary = native_roi, rsc = read_csv(baseline$rsc_csv),
    fate = read_csv(baseline$fate_csv),
    endpoint_summary = read_csv(baseline$endpoint_summary_csv),
    tractogram_density_summary = contexts[[1L]]$summary,
    report_indexes = list(native_tck = native$report_index, equivalent_trk = equivalent$report_index),
    parity_passed = all(parity$passed), validation_summary = summarize_validation(),
    tolerances = state$tolerances,
    interpretation = state$scope
  )
  jsonlite::write_json(science, file.path(output_dir, "science-summary.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  writeLines(c(
    "Native TCK public-data validation", "",
    paste("Native TCK reports:", native$report_index),
    paste("Equivalent TRK reports:", equivalent$report_index), "",
    "validation-results.json records software, input hashes, settings, timings, warnings and checks.",
    "parity-checks.csv records exact comparisons and the continuous-coordinate tolerance.",
    "Native TCK execution/invariance and TCK-versus-TRK equivalence have separate status summaries.",
    "Overall validation fails if any strict comparison fails; separate summaries do not waive differences.",
    "science-summary.json contains manuscript-traceable results and interpretation limits.",
    "The density-context folders contain independently generated streamline-occupancy maps.",
    "Four analysis runs use one dataset and one subject; these are computational repeats.",
    "The masks are tissue/RAS-x geometric demonstrations, not validated atlas or hemisphere regions.",
    "A computational pass does not replace visual anatomical review."
  ), file.path(output_dir, "review-files.txt"))
  if (!all(parity$passed)) {
    stop("Parity failed; inspect parity-checks.csv. No differing results were adjusted.", call. = FALSE)
  }
}

write_state()
save_session <- function() writeLines(capture.output(utils::sessionInfo()),
                                      file.path(output_dir, "session-info.txt"))
tryCatch(
  withCallingHandlers({
    main()
    state$status <- "passed"
    state$completed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    write_state()
    save_session()
    message("Validation passed: ", output_dir)
  }, warning = function(w) {
    state$warnings <<- c(state$warnings, conditionMessage(w))
    write_state()
  }),
  error = function(e) {
    state$status <<- "failed"
    state$error <<- conditionMessage(e)
    state$completed_at <<- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    write_state()
    save_session()
    stop(conditionMessage(e), call. = FALSE)
  }
)
