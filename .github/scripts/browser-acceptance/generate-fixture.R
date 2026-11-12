#!/usr/bin/env Rscript
# Use an installed candidate and only the package-owned synthetic fixture.
# Usage: Rscript generate-fixture.R <new-output-directory> [trk|tck]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L, 2L))
format <- if (length(args) == 2L) match.arg(args[[2L]], c("trk", "tck")) else "trk"
root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
if (dir.exists(root)) stop("Choose a new browser acceptance output directory.")
dir.create(root, recursive = TRUE)
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
original <- file.path(root, "original report # %20 caf\u00e9")
relocated <- file.path(root, "relocated report # %20 caf\u00e9")

library(FibrTraQ)
run_tck_example <- function() {
  data <- ftq_example_data(verify = TRUE)
  params <- data$parameters
  stopifnot(identical(params$track_space, "voxel"),
            identical(as.logical(params$axis_flip), rep(FALSE, 3L)))
  input_dir <- file.path(root, "native-tck-inputs")
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
  tracks <- freesurferformats::read.dti.trk(data$tractogram)$tracks
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
      writeBin(as.numeric(t(track$coords - 1)), con, size = 4L, endian = "little")
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

example <- if (format == "trk") {
  ftq_run_example(output_dir = original, backend = "cpp", verbose = FALSE)
} else {
  run_tck_example()
}
stopifnot(isTRUE(example$validation$passed),
          nrow(example$validation$tables) == 7L,
          all(example$validation$tables$passed),
          nrow(example$validation$claims) == 18L,
          all(example$validation$claims$passed),
          all(example$validation$reports$passed))
original <- normalizePath(original, winslash = "/", mustWork = TRUE)
files <- list.files(original, recursive = TRUE, all.files = TRUE,
                    no.. = TRUE, full.names = FALSE)
files <- files[!dir.exists(file.path(original, files))]
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
