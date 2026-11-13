#!/usr/bin/env Rscript
# Use an installed candidate and only the package-owned synthetic fixture.
# Usage: Rscript generate-fixture.R <new-output-directory> [trk|tck] [normal|near-limit]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 1:3)
format <- if (length(args) >= 2L) match.arg(args[[2L]], c("trk", "tck")) else "trk"
path_mode <- if (length(args) == 3L) match.arg(args[[3L]], c("normal", "near-limit")) else "normal"
root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
if (dir.exists(root)) stop("Choose a new browser acceptance output directory.")
dir.create(root, recursive = TRUE)
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
original <- file.path(root, "original report # %20 caf\u00e9")
relocated <- file.path(root, "relocated report # %20 caf\u00e9")

library(FibrTraQ)
run_tck_example <- function(original, input_tag) {
  data <- ftq_example_data(verify = TRUE)
  params <- data$parameters
  stopifnot(identical(params$track_space, "voxel"),
            identical(as.logical(params$axis_flip), rep(FALSE, 3L)))
  input_dir <- file.path(root, paste0("native-tck-inputs-", input_tag))
  dir.create(input_dir)

  # This package-owned fixture declares one-based voxel coordinates. Give its
  # copied images an explicit identity affine in millimetres; then each voxel
  # centre (i,j,k) has TCK world coordinates (i-1,j-1,k-1). The bundled inputs
  # remain immutable, and no package coordinate-mapping function writes TCK.
  image_copy <- function(path) {
    values <- as.array(RNifti::readNifti(path))
    image <- RNifti::asNifti(values, reference = list(
      sform_code = 2L, qform_code = 0L,
      pixdim = c(1, 1, 1, 1, 1, 1, 1, 1),
      srow_x = c(1, 0, 0, 0), srow_y = c(0, 1, 0, 0),
      srow_z = c(0, 0, 1, 0)
    ))
    RNifti::pixunits(image) <- "mm"
    destination <- file.path(input_dir, basename(path))
    RNifti::writeNifti(image, destination)
    destination
  }
  focal <- image_copy(data$focal_mask)
  targets <- vapply(data$target_masks, image_copy, character(1L))
  local <- image_copy(data$local_mask)
  reference <- image_copy(data$reference_image)
  contexts <- vapply(data$endpoint_contexts, image_copy, character(1L))
  read_trk <- freesurferformats::read.dti.trk
  read_args <- list(filepath = data$tractogram)
  if ("coords" %in% names(formals(read_trk))) read_args$coords <- "native"
  tracks <- as.list(do.call(read_trk, read_args)$tracks)
  tractogram <- file.path(input_dir, "sub-example_tractogram.tck")
  write_tck <- function() {
    offset <- 512L
    header <- charToRaw(paste0(
      "mrtrix tracks\ncount: ", length(tracks),
      "\ndatatype: Float32LE\nfile: . ", offset, "\nEND\n"
    ))
    stopifnot(length(header) < offset)
    con <- file(tractogram, "wb")
    on.exit(close(con), add = TRUE)
    writeBin(c(header, raw(offset - length(header))), con)
    for (track in tracks) {
      # freesurferformats versions return either coordinate matrices or
      # list records. Normalize the representation without changing points.
      coords <- as_track(track)$coords
      writeBin(as.numeric(t(coords - 1)), con, size = 4L, endian = "little")
      writeBin(rep(NaN, 3L), con, size = 4L, endian = "little")
    }
    writeBin(rep(Inf, 3L), con, size = 4L, endian = "little")
  }
  write_tck()
  analysis <- ftq_run_files(
    tractogram = tractogram, focal_mask = focal, target_masks = targets,
    output_dir = original, subject = params$subject, focal_label = params$focal_label,
    local_mask = local, reference_image = reference, endpoint_contexts = contexts,
    description = params$description, track_space = "auto",
    discretization = params$discretization,
    focal_assignment_mode = params$focal_assignment_mode,
    streamline_voxelization_mode = params$streamline_voxelization_mode,
    assignment_mode = params$target_intersection_mode,
    endpoint_analysis = params$endpoint_analysis, track_reader = params$track_reader,
    track_chunk_size = as.integer(params$track_chunk_size), backend = "cpp",
    write_maps = TRUE, write_report = TRUE, write_exploratory_report = TRUE,
    write_qc_viewer = TRUE,
    sample_per_target = as.integer(params$example_execution$sample_per_target),
    verbose = FALSE, progress = FALSE, open = FALSE
  )
  run <- analysis$index[1L, , drop = FALSE]
  actual <- stats::setNames(vapply(names(data$expected), function(artifact) {
    as.character(run[[paste0(artifact, "_csv")]])
  }, character(1L)), names(data$expected))

  # Only coordinate provenance and the copied context-image checksum differ.
  # Keep every scientific value, column and row, including all 18 claims.
  expected <- data$expected
  roi <- utils::read.csv(expected[["roi_summary"]], check.names = FALSE)
  stopifnot(identical(roi$selected_track_space, "voxel"),
            identical(roi$track_voxel_size_mm, "1,1,1"),
            identical(roi$mapping_selection_rule, "pinned"))
  roi$selected_track_space <- "world"
  roi$track_voxel_size_mm <- NA_character_
  roi$mapping_selection_rule <- "tck_world_ras_mm_from_nifti_affine"
  expected[["roi_summary"]] <- file.path(input_dir, "expected-tck-roi-summary.csv")
  utils::write.csv(roi, expected[["roi_summary"]], row.names = FALSE)
  endpoint_context <- utils::read.csv(expected[["endpoint_context"]], check.names = FALSE)
  context_rows <- match(names(contexts), endpoint_context$context)
  stopifnot(!anyNA(context_rows))
  endpoint_context$checksum_md5[context_rows] <- unname(tools::md5sum(contexts))
  expected[["endpoint_context"]] <- file.path(input_dir, "expected-tck-endpoint-context.csv")
  utils::write.csv(endpoint_context, expected[["endpoint_context"]], row.names = FALSE)
  tables <- getFromNamespace("ftq_example_compare_tables", "FibrTraQ")(actual, expected)
  claims <- getFromNamespace("ftq_example_check_claims", "FibrTraQ")(
    actual, expected, data$claims
  )
  report_paths <- c(
    report_index = analysis$report_index, primary_report = run$primary_report,
    exploratory_report = run$exploratory_report, interactive_3d_qc = run$interactive_3d_qc
  )
  reports <- data.frame(
    artifact = names(report_paths), path = unname(report_paths),
    passed = file.exists(report_paths) & file.info(report_paths)$size > 0,
    details = "generated and non-empty", stringsAsFactors = FALSE
  )
  validation <- list(
    passed = isTRUE(run$validation_passed) && all(tables$passed) &&
      all(claims$passed) && all(reports$passed),
    package_validation_passed = isTRUE(run$validation_passed),
    tables = tables, claims = claims, reports = reports
  )
  getFromNamespace("ftq_write_regional_study_index", "FibrTraQ")(
    run_index = analysis$index, path = analysis$report_index,
    output_dir = original, description = params$description,
    example_validation = validation
  )
  # Recheck immutable source files after generating the independent encoding.
  ftq_example_data(verify = TRUE)
  c(list(analysis = analysis, validation = validation), as.list(report_paths))
}

run_example <- function(output_dir, input_tag) {
  if (format == "trk") {
    ftq_run_example(output_dir = output_dir, backend = "cpp", verbose = FALSE)
  } else {
    run_tck_example(output_dir, input_tag)
  }
}
check_example <- function(example) {
  stopifnot(isTRUE(example$validation$passed),
          nrow(example$validation$tables) == 7L,
          all(example$validation$tables$passed),
          nrow(example$validation$claims) == 18L,
          all(example$validation$claims$passed),
          all(example$validation$reports$passed))
}
report_files <- function(directory) {
  files <- list.files(directory, recursive = TRUE, all.files = TRUE,
                      no.. = TRUE, full.names = FALSE)
  files[!dir.exists(file.path(directory, files))]
}
utf16_length <- function(paths) {
  vapply(iconv(enc2utf8(paths), from = "UTF-8", to = "UTF-16LE", toRaw = TRUE),
         function(value) {
           if (is.null(value)) stop("Could not measure a fixture path in UTF-16.")
           length(value) / 2
         }, numeric(1L), USE.NAMES = FALSE)
}
planned_paths <- function(output_dir) {
  data <- ftq_example_data(verify = TRUE)
  params <- data$parameters
  getFromNamespace("ftq_regional_study_output_paths", "FibrTraQ")(
    output_dir, subjects = params$subject,
    region_labels = c(params$focal_label, names(data$target_masks)),
    focal_regions = params$focal_label, description = params$description,
    write_endpoint = identical(params$endpoint_analysis, "exact_and_dilated1"),
    write_local_search = identical(params$target_intersection_mode, "local_search")
  )
}
example <- run_example(original, "baseline")
check_example(example)
original <- normalizePath(original, winslash = "/", mustWork = TRUE)
files <- report_files(original)
path_boundary <- NULL
if (path_mode == "near-limit") {
  stopifnot(.Platform$OS.type == "windows")
  baseline_root <- original
  baseline_files <- files
  baseline_checksums <- unname(tools::md5sum(file.path(original, files)))
  longest_relative <- max(utf16_length(files))
  sized_root <- function(prefix, maximum_file_length) {
    desired_root_length <- maximum_file_length - longest_relative - 1L
    padding <- desired_root_length - utf16_length(prefix)
    if (padding < 1L) stop("Choose a shorter browser acceptance root for near-limit checks.")
    paste0(prefix, strrep("p", padding))
  }
  # Run the writers again at the measured boundary; a relocated short-root run
  # alone cannot demonstrate that R creates report images and tables here.
  original <- sized_root(file.path(root, "cran-temp", "generated # %20 caf\u00e9 "), 255L)
  relocated <- sized_root(file.path(root, "relocated # %20 caf\u00e9 "), 259L)
  original_planned <- planned_paths(original)
  example <- run_example(original, "near-limit")
  check_example(example)
  original <- normalizePath(original, winslash = "/", mustWork = TRUE)
  files <- report_files(original)
  stopifnot(identical(files, baseline_files),
            max(utf16_length(file.path(original, files))) == 255L,
            identical(unname(tools::md5sum(file.path(baseline_root, baseline_files))),
                      baseline_checksums))

  rejected_root <- sized_root(file.path(root, "cran-temp", "rejected "), 300L)
  rejected_planned <- planned_paths(rejected_root)
  sentinel <- file.path(root, "existing-output-sentinel.txt")
  writeLines("Existing output must remain intact after a rejected request.", sentinel)
  sentinel_checksum <- unname(tools::md5sum(sentinel))
  rejection <- tryCatch({ run_example(rejected_root, "over-limit"); NULL }, error = identity)
  # Write the observed outcome before assertions so a broken guard still leaves
  # enough evidence to distinguish a rejection from a partial-output failure.
  rejected_evidence <- list(
    output_root = rejected_root,
    planned_maximum_utf16_units = max(utf16_length(rejected_planned)),
    message = if (inherits(rejection, "error")) conditionMessage(rejection) else NULL,
    output_root_absent = !file.exists(rejected_root),
    existing_output_sentinel = sentinel,
    existing_output_sentinel_md5 = sentinel_checksum,
    existing_output_sentinel_preserved = identical(unname(tools::md5sum(sentinel)), sentinel_checksum),
    baseline_outputs_preserved = identical(
      unname(tools::md5sum(file.path(baseline_root, baseline_files))), baseline_checksums)
  )
  jsonlite::write_json(rejected_evidence, file.path(root, "over-limit-evidence.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  stopifnot(inherits(rejection, "error"),
            startsWith(conditionMessage(rejection),
              "Windows output path exceeds the supported limit of 259 UTF-16 code units"),
            rejected_evidence$output_root_absent,
            rejected_evidence$existing_output_sentinel_preserved,
            rejected_evidence$baseline_outputs_preserved)
  path_boundary <- list(
    baseline_root = baseline_root,
    generated_directly_at_boundary = TRUE,
    original_expected_maximum_utf16_units = 255L,
    original_planned_maximum_utf16_units = max(utf16_length(original_planned)),
    relocated_expected_maximum_utf16_units = 259L,
    relative_inventory_identical = TRUE,
    rejected_request = rejected_evidence
  )
}
dir.create(relocated)
for (relative in files) {
  destination <- file.path(relocated, relative)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(original, relative), destination, overwrite = FALSE)) {
    stop("Could not copy report asset: ", relative)
  }
}
stopifnot(identical(unname(tools::md5sum(file.path(original, files))),
                    unname(tools::md5sum(file.path(relocated, files)))))
relocated <- normalizePath(relocated, winslash = "/", mustWork = TRUE)
if (path_mode == "near-limit") {
  stopifnot(max(utf16_length(file.path(relocated, files))) == 259L)
}
relative_path <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  stopifnot(startsWith(path, paste0(original, "/")))
  substring(path, nchar(original) + 2L)
}
reports <- list(
  index = relative_path(example$report_index),
  results = relative_path(example$primary_report),
  target = relative_path(example$exploratory_report),
  viewer = relative_path(example$interactive_3d_qc)
)
metadata <- example$analysis$index$metadata_json[[1L]]
manifest <- list(
  fixture = paste0("package-owned synthetic known-answer example (", toupper(format),
                   "); not anatomical MRI"),
  tractogram_format = format,
  path_mode = path_mode,
  path_boundary = path_boundary,
  tck_coordinate_convention = if (format == "tck") "RAS+ world millimetres" else NULL,
  tck_image_geometry = if (format == "tck") "identity sform; millimetre units" else NULL,
  package_version = as.character(packageVersion("FibrTraQ")),
  source_commit = Sys.getenv("GITHUB_SHA", unset = "local-uncommitted-review"),
  r_version = R.version.string,
  os = unname(Sys.info()[["sysname"]]),
  original_root = original,
  relocated_root = relocated,
  reports = reports,
  metadata = relative_path(metadata),
  html_files = as.list(files[grepl("[.]html$", files, ignore.case = TRUE)]),
  file_count = length(files),
  copy_checksums_identical = TRUE,
  scientific_tables_passed = 7L,
  declared_values_passed = 18L,
  original_outputs_preserved = TRUE
)
jsonlite::write_json(manifest, file.path(root, "browser-fixture.json"),
                     auto_unbox = TRUE, pretty = TRUE)
cat("Browser fixture generated with ", length(files), " files.\n", sep = "")
