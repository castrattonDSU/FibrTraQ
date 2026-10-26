# Regenerate the package-owned FibrTraQ known-answer example.
#
# Run from the repository root after loading the current package source:
#   Rscript data-raw/generate-known-answer-example.R
#
# This developer script creates immutable example inputs, runs the reference R
# backend, verifies independently specified counts, and snapshots only the
# scientific CSV outputs. Reports and temporary run metadata are not bundled.

if (!file.exists("DESCRIPTION") ||
    !grepl("^Package: FibrTraQ$", readLines("DESCRIPTION", n = 1L))) {
  stop("Run this script from the FibrTraQ repository root.", call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE) ||
    !requireNamespace("RNifti", quietly = TRUE) ||
    !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("devtools, RNifti, and jsonlite are required to regenerate the example.",
       call. = FALSE)
}

devtools::load_all(".", quiet = TRUE)

fixture_root <- file.path("inst", "extdata", "fibrtraq-known-answer")
input_dir <- file.path(fixture_root, "inputs")
expected_dir <- file.path(fixture_root, "expected")
unlink(fixture_root, recursive = TRUE, force = TRUE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(expected_dir, recursive = TRUE, showWarnings = FALSE)

write_fixed <- function(con, text, n) {
  bytes <- rep(as.raw(0), n)
  text_bytes <- charToRaw(text)
  if (length(text_bytes)) {
    keep <- seq_len(min(length(text_bytes), n))
    bytes[keep] <- text_bytes[keep]
  }
  writeBin(bytes, con)
}

canonicalize_lf <- function(path) {
  bytes <- readBin(
    path,
    what = "raw",
    n = as.integer(file.info(path)$size)
  )
  text <- rawToChar(bytes)
  text <- gsub("\r\n", "\n", text, fixed = TRUE)
  text <- gsub("\r", "\n", text, fixed = TRUE)
  writeBin(charToRaw(text), path)
  invisible(path)
}

write_trk <- function(path, tracks, dims) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  endian <- "little"

  write_fixed(con, "TRACK", 6L)
  writeBin(as.integer(dims), con, size = 2L, endian = endian)
  writeBin(as.numeric(c(1, 1, 1)), con, size = 4L, endian = endian)
  writeBin(as.numeric(c(0, 0, 0)), con, size = 4L, endian = endian)
  writeBin(as.integer(0), con, size = 2L, endian = endian)
  write_fixed(con, "", 200L)
  writeBin(as.integer(0), con, size = 2L, endian = endian)
  write_fixed(con, "", 200L)
  writeBin(as.numeric(t(diag(4))), con, size = 4L, endian = endian)
  write_fixed(con, "", 444L)
  write_fixed(con, "RAS", 4L)
  write_fixed(con, "", 4L)
  writeBin(as.numeric(c(1, 0, 0, 0, 1, 0)), con,
           size = 4L, endian = endian)
  write_fixed(con, "", 2L)
  writeBin(as.integer(rep(0, 6)), con, size = 1L, endian = endian)
  writeBin(as.integer(length(tracks)), con, size = 4L, endian = endian)
  writeBin(as.integer(2), con, size = 4L, endian = endian)
  writeBin(as.integer(1000), con, size = 4L, endian = endian)

  for (track in tracks) {
    coords <- as.matrix(track)[, 1:3, drop = FALSE]
    writeBin(as.integer(nrow(coords)), con, size = 4L, endian = endian)
    writeBin(as.numeric(t(coords)), con, size = 4L, endian = endian)
  }
  invisible(path)
}

dims <- c(9L, 9L, 9L)
blank <- array(0, dim = dims)
focal <- blank
cortex <- blank
brainstem <- blank
local <- blank
context <- array(1, dim = dims)

focal[2, 4, 4] <- 1
focal[2, 5, 4] <- 1
cortex[6:8, 4, 4] <- 1
brainstem[2, 7:8, 4] <- 1
local[1:3, 3:6, 3:5] <- 1

reference <- array(0, dim = dims)
for (z in seq_len(dims[3])) {
  for (y in seq_len(dims[2])) {
    for (x in seq_len(dims[1])) {
      reference[x, y, z] <- 100 + (10 * x) + (3 * y) + z
    }
  }
}

paths <- list(
  tractogram = file.path(input_dir, "sub-example_tractogram.trk"),
  focal = file.path(input_dir, "sub-example_focal.nii.gz"),
  cortex = file.path(input_dir, "sub-example_target-cortex.nii.gz"),
  brainstem = file.path(input_dir, "sub-example_target-brainstem.nii.gz"),
  local = file.path(input_dir, "sub-example_local.nii.gz"),
  context = file.path(input_dir, "sub-example_context.nii.gz"),
  reference = file.path(input_dir, "sub-example_reference.nii.gz")
)

RNifti::writeNifti(focal, paths$focal)
template <- RNifti::readNifti(paths$focal, internal = TRUE)
RNifti::writeNifti(cortex, paths$cortex, template = template)
RNifti::writeNifti(brainstem, paths$brainstem, template = template)
RNifti::writeNifti(local, paths$local, template = template)
RNifti::writeNifti(context, paths$context, template = template)
RNifti::writeNifti(reference, paths$reference, template = template)

tracks <- list(
  # Cortex: two exact focal-target endpoint pairs.
  rbind(
    c(2, 4, 4), c(3, 4, 4), c(4, 4, 4), c(5, 4, 4),
    c(6, 4, 4), c(7, 4, 4), c(8, 4, 4)
  ),
  rbind(
    c(8, 4, 4), c(7, 4, 4), c(6, 4, 4), c(5, 4, 4),
    c(4, 4, 4), c(3, 4, 4), c(2, 4, 4)
  ),
  # Cortex: focal-only, target-only, neither, and dilated-pair classes.
  rbind(
    c(2, 5, 4), c(3, 5, 4), c(4, 5, 4), c(5, 5, 4),
    c(6, 4, 4), c(7, 4, 4), c(8, 6, 4), c(9, 7, 4)
  ),
  rbind(
    c(1, 1, 4), c(1, 2, 4), c(1, 3, 4), c(2, 4, 4),
    c(3, 4, 4), c(4, 4, 4), c(5, 4, 4), c(6, 4, 4)
  ),
  rbind(
    c(1, 1, 4), c(1, 2, 4), c(2, 5, 4), c(3, 5, 4),
    c(4, 5, 4), c(5, 5, 4), c(6, 4, 4), c(7, 4, 4),
    c(8, 6, 4), c(9, 7, 4)
  ),
  rbind(
    c(1, 4, 4), c(2, 4, 4), c(3, 4, 4), c(4, 4, 4),
    c(5, 4, 4), c(6, 4, 4), c(7, 4, 4), c(8, 5, 4)
  ),
  # Brainstem: two exact pairs, focal-only, and target-only classes.
  rbind(c(2, 4, 4), c(2, 5, 4), c(2, 6, 4), c(2, 7, 4), c(2, 8, 4)),
  rbind(c(2, 5, 4), c(2, 6, 4), c(2, 7, 4), c(2, 8, 4)),
  rbind(
    c(2, 4, 4), c(2, 5, 4), c(2, 6, 4), c(2, 7, 4),
    c(3, 8, 4), c(4, 9, 4)
  ),
  rbind(
    c(1, 1, 4), c(1, 2, 4), c(1, 3, 4), c(2, 4, 4),
    c(2, 5, 4), c(2, 6, 4), c(2, 7, 4), c(2, 8, 4)
  ),
  # Focal streamlines with no target, including a one-point streamline.
  rbind(c(2, 4, 4), c(3, 3, 4), c(4, 2, 4)),
  matrix(c(2, 5, 4), nrow = 1L),
  # Non-focal controls.
  rbind(c(6, 4, 4), c(7, 4, 4), c(8, 4, 4)),
  rbind(c(1, 1, 1), c(1, 2, 1))
)
write_trk(paths$tractogram, tracks, dims)

parameters <- list(
  dataset_id = "fibrtraq-known-answer-v1",
  dataset_type = "package-authored synthetic known-answer fixture",
  biological_interpretation = FALSE,
  subject = "sub-example",
  focal_label = "focal",
  target_labels = c("cortex", "brainstem"),
  description = "known-answer",
  track_space = "voxel",
  discretization = "round",
  axis_flip = c(FALSE, FALSE, FALSE),
  focal_assignment_mode = "all_focal_voxels",
  streamline_voxelization_mode = "sample_points",
  target_intersection_mode = "full_streamline",
  endpoint_analysis = "exact_and_dilated1",
  endpoint_context_labels = "analysis_domain",
  track_reader = "stream",
  track_chunk_size = 4L,
  expected_table_backend = "r",
  example_execution = list(
    write_maps = TRUE,
    write_report = TRUE,
    write_exploratory_report = TRUE,
    write_qc_viewer = TRUE,
    sample_per_target = 6L,
    progress = FALSE
  ),
  expected_table_generation = list(
    write_maps = FALSE,
    write_report = FALSE,
    write_exploratory_report = FALSE,
    write_qc_viewer = FALSE,
    progress = FALSE
  ),
  generation_provenance = list(
    generator = "data-raw/generate-known-answer-example.R",
    package_version = unname(read.dcf("DESCRIPTION")[1L, "Version"]),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    RNifti_version = as.character(utils::packageVersion("RNifti")),
    jsonlite_version = as.character(utils::packageVersion("jsonlite")),
    devtools_version = as.character(utils::packageVersion("devtools"))
  )
)
jsonlite::write_json(
  parameters,
  file.path(fixture_root, "analysis-parameters.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
writeLines(
  c(
    "# FibrTraQ Known-Answer Example",
    "",
    "This package-authored synthetic dataset provides a small, immutable",
    "end-to-end example for FibrTraQ. It contains 14 TrackVis streamlines on",
    "a 9 x 9 x 9 identity-affine voxel grid, a two-voxel focal mask, cortex",
    "and brainstem target masks, a local mask, an endpoint context mask, and",
    "a reference image.",
    "",
    "The dataset is designed to exercise focal-voxel assignment, full-path",
    "target intersections, mutually exclusive fate categories, exact and",
    "one-voxel endpoint classes, target penetration, and report generation.",
    "It is not anatomical MRI and must not be used as biological evidence.",
    "",
    "The `expected` directory contains seven scientific CSV tables generated",
    "with the reference R backend after independent count assertions in the",
    "[fixture generator](https://github.com/castrattonDSU/FibrTraQ/blob/main/data-raw/generate-known-answer-example.R).",
    "`ftq_run_example()` compares new outputs against those tables and resolves",
    "all 18 declarations in `known-answer-claims.csv` from both the stored and",
    "newly generated tables. Backend identity and installation-path provenance",
    "are normalized before comparison; all scientific values, labels,",
    "denominators, and row ordering must match exactly.",
    "",
    "The fixture was authored for FibrTraQ and is distributed under the same",
    "GNU General Public License, version 2 or later, as the package."
  ),
  file.path(fixture_root, "README.md")
)

run_root <- tempfile("fibrtraq-known-answer-generation-")
on.exit(unlink(run_root, recursive = TRUE, force = TRUE), add = TRUE)
analysis <- ftq_run_files(
  tractogram = paths$tractogram,
  focal_mask = paths$focal,
  target_masks = c(cortex = paths$cortex, brainstem = paths$brainstem),
  output_dir = run_root,
  subject = parameters$subject,
  focal_label = parameters$focal_label,
  local_mask = paths$local,
  reference_image = paths$reference,
  endpoint_contexts = c(analysis_domain = paths$context),
  description = parameters$description,
  track_space = parameters$track_space,
  discretization = parameters$discretization,
  axis_flip = parameters$axis_flip,
  focal_assignment_mode = parameters$focal_assignment_mode,
  streamline_voxelization_mode = parameters$streamline_voxelization_mode,
  endpoint_analysis = parameters$endpoint_analysis,
  track_reader = parameters$track_reader,
  track_chunk_size = parameters$track_chunk_size,
  backend = parameters$expected_table_backend,
  write_maps = parameters$expected_table_generation$write_maps,
  write_report = parameters$expected_table_generation$write_report,
  write_exploratory_report =
    parameters$expected_table_generation$write_exploratory_report,
  write_qc_viewer = parameters$expected_table_generation$write_qc_viewer,
  overwrite = TRUE,
  verbose = FALSE,
  progress = parameters$expected_table_generation$progress
)

run <- analysis$index[1, , drop = FALSE]
result <- utils::read.csv(run$result_csv, check.names = FALSE)
rsc <- utils::read.csv(run$rsc_csv, check.names = FALSE)
fate <- utils::read.csv(run$fate_csv, check.names = FALSE)
roi <- utils::read.csv(run$roi_summary_csv, check.names = FALSE)
endpoint <- utils::read.csv(run$endpoint_summary_csv, check.names = FALSE)
penetration <- utils::read.csv(
  run$endpoint_penetration_csv,
  check.names = FALSE
)
context_summary <- utils::read.csv(
  run$endpoint_context_csv,
  check.names = FALSE
)

stopifnot(
  identical(result$fNode_Coords, c("2_4_4", "2_5_4")),
  identical(as.integer(result$tracks), c(8L, 7L)),
  identical(as.integer(result$local), c(0L, 1L)),
  identical(as.integer(result$leave_no_target), c(1L, 0L)),
  identical(as.integer(result$cortex), c(4L, 2L)),
  identical(as.integer(result$brainstem), c(3L, 4L)),
  identical(as.integer(roi$processed_streamlines), 14L),
  identical(as.integer(roi$unique_focal_intersecting_streamlines), 12L),
  identical(as.integer(roi$unique_local_streamlines), 1L),
  identical(as.integer(roi$unique_leave_no_target_streamlines), 1L),
  identical(as.integer(roi$unique_leave_with_target_streamlines), 10L),
  identical(as.integer(roi$total_focal_voxel_streamline_associations), 15L),
  identical(as.integer(roi[["unique_target_streamlines:cortex"]]), 6L),
  identical(as.integer(roi[["unique_target_streamlines:brainstem"]]), 4L),
  identical(rsc$node, c("cortex", "brainstem")),
  identical(as.integer(rsc$target_associations), c(6L, 7L)),
  identical(as.integer(rsc$voxels_connected), c(2L, 2L)),
  identical(as.numeric(rsc$rsc_percent), c(100, 100)),
  identical(fate$fate, c("local", "leave_no_target", "leave_with_target")),
  identical(as.integer(fate$tracks), c(1L, 1L, 13L))
)

cortex_endpoint <- endpoint[endpoint$target == "cortex", , drop = FALSE]
cortex_exact <- cortex_endpoint[cortex_endpoint$tolerance == "exact", ]
cortex_dilated <- cortex_endpoint[cortex_endpoint$tolerance == "dilated1", ]
brainstem_exact <- endpoint[
  endpoint$target == "brainstem" & endpoint$tolerance == "exact",
]
brainstem_dilated <- endpoint[
  endpoint$target == "brainstem" & endpoint$tolerance == "dilated1",
]
cortex_penetration <- penetration[
  penetration$target == "cortex",
  ,
  drop = FALSE
]
brainstem_penetration <- penetration[
  penetration$target == "brainstem",
  ,
  drop = FALSE
]
grid_context <- context_summary[
  context_summary$context == "reference_grid",
  ,
  drop = FALSE
]
mask_context <- context_summary[
  context_summary$context == "analysis_domain",
  ,
  drop = FALSE
]
stopifnot(
  nrow(endpoint) == 4L,
  cortex_exact$target_intersection_streamlines == 6,
  cortex_exact$focal_endpoint_streamlines == 8,
  cortex_exact$paired_focal_target_endpoints == 2,
  cortex_exact$target_endpoint_only == 1,
  cortex_exact$focal_endpoint_only == 1,
  cortex_exact$neither_relevant_endpoint == 2,
  cortex_exact$shared_mask_ambiguous == 0,
  cortex_dilated$focal_endpoint_streamlines == 9,
  cortex_dilated$paired_focal_target_endpoints == 3,
  cortex_dilated$target_endpoint_only == 1,
  cortex_dilated$focal_endpoint_only == 1,
  cortex_dilated$neither_relevant_endpoint == 1,
  cortex_dilated$shared_mask_ambiguous == 0,
  brainstem_exact$target_intersection_streamlines == 4,
  brainstem_exact$focal_endpoint_streamlines == 8,
  brainstem_exact$paired_focal_target_endpoints == 2,
  brainstem_exact$target_endpoint_only == 1,
  brainstem_exact$focal_endpoint_only == 1,
  brainstem_exact$neither_relevant_endpoint == 0,
  brainstem_dilated$focal_endpoint_streamlines == 9,
  brainstem_dilated$paired_focal_target_endpoints == 2,
  identical(
    as.integer(cortex_penetration$unique_target_voxels),
    c(1L, 2L, 3L)
  ),
  identical(
    as.integer(cortex_penetration$streamline_count),
    c(1L, 3L, 2L)
  ),
  identical(
    as.integer(brainstem_penetration$unique_target_voxels),
    c(1L, 2L)
  ),
  identical(
    as.integer(brainstem_penetration$streamline_count),
    c(1L, 3L)
  ),
  grid_context$valid_streamlines == 14,
  grid_context$endpoint_count == 27,
  grid_context$endpoints_inside == 27,
  grid_context$endpoints_outside == 0,
  grid_context$streamlines_with_both_endpoints_inside == 14,
  mask_context$endpoint_count == 27,
  mask_context$endpoints_inside == 27,
  mask_context$endpoints_outside == 0
)

expected_sources <- c(
  result = run$result_csv,
  rsc = run$rsc_csv,
  fate = run$fate_csv,
  roi_summary = run$roi_summary_csv,
  endpoint_summary = run$endpoint_summary_csv,
  endpoint_penetration = run$endpoint_penetration_csv,
  endpoint_context = run$endpoint_context_csv
)
expected_names <- paste0(names(expected_sources), ".csv")
copied <- file.copy(
  unname(expected_sources),
  file.path(expected_dir, expected_names),
  overwrite = TRUE,
  copy.mode = FALSE
)
if (!all(copied)) {
  stop("Failed to copy one or more expected scientific tables.", call. = FALSE)
}
expected_context_path <- file.path(expected_dir, "endpoint_context.csv")
expected_context <- utils::read.csv(
  expected_context_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
has_source <- !is.na(expected_context$source_path) &
  nzchar(expected_context$source_path)
expected_context$source_path[has_source] <- basename(
  expected_context$source_path[has_source]
)
utils::write.csv(expected_context, expected_context_path, row.names = FALSE)

claims <- data.frame(
  claim = c(
    "processed_streamlines",
    "unique_focal_intersecting_streamlines",
    "total_focal_voxel_streamline_associations",
    "unique_any_target_positive_streamlines",
    "unique_local_streamlines",
    "unique_leave_no_target_streamlines",
    "cortex_target_intersections",
    "brainstem_target_intersections",
    "cortex_focal_voxel_target_associations",
    "brainstem_focal_voxel_target_associations",
    "cortex_exact_paired_endpoints",
    "cortex_dilated1_paired_endpoints",
    "brainstem_exact_paired_endpoints",
    "reference_grid_endpoint_count",
    "reference_grid_endpoints_inside",
    "fate_local_associations",
    "fate_leave_no_target_associations",
    "fate_leave_with_target_associations"
  ),
  source_table = c(
    rep("roi_summary", 6),
    rep("endpoint_summary", 2),
    rep("rsc", 2),
    rep("endpoint_summary", 3),
    rep("endpoint_context", 2),
    rep("fate", 3)
  ),
  value_column = c(
    "processed_streamlines",
    "unique_focal_intersecting_streamlines",
    "total_focal_voxel_streamline_associations",
    "unique_any_target_positive_streamlines",
    "unique_local_streamlines",
    "unique_leave_no_target_streamlines",
    "target_intersection_streamlines",
    "target_intersection_streamlines",
    "target_associations",
    "target_associations",
    "paired_focal_target_endpoints",
    "paired_focal_target_endpoints",
    "paired_focal_target_endpoints",
    "endpoint_count",
    "endpoints_inside",
    "tracks",
    "tracks",
    "tracks"
  ),
  filter_column_1 = c(
    rep(NA_character_, 6),
    "target",
    "target",
    "node",
    "node",
    "target",
    "target",
    "target",
    "context",
    "context",
    "fate",
    "fate",
    "fate"
  ),
  filter_value_1 = c(
    rep(NA_character_, 6),
    "cortex",
    "brainstem",
    "cortex",
    "brainstem",
    "cortex",
    "cortex",
    "brainstem",
    "reference_grid",
    "reference_grid",
    "local",
    "leave_no_target",
    "leave_with_target"
  ),
  filter_column_2 = c(
    rep(NA_character_, 6),
    "tolerance",
    "tolerance",
    rep(NA_character_, 2),
    "tolerance",
    "tolerance",
    "tolerance",
    rep(NA_character_, 5)
  ),
  filter_value_2 = c(
    rep(NA_character_, 6),
    "exact",
    "exact",
    rep(NA_character_, 2),
    "exact",
    "dilated1",
    "exact",
    rep(NA_character_, 5)
  ),
  expected_value = c(
    14, 12, 15, 10, 1, 1, 6, 4, 6, 7, 2, 3, 2, 27, 27, 1, 1, 13
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  claims,
  file.path(fixture_root, "known-answer-claims.csv"),
  row.names = FALSE
)

text_paths <- c(
  file.path(fixture_root, "analysis-parameters.json"),
  file.path(fixture_root, "known-answer-claims.csv"),
  file.path(fixture_root, "README.md"),
  list.files(expected_dir, pattern = "\\.csv$", full.names = TRUE)
)
invisible(lapply(text_paths, canonicalize_lf))

inventory_paths <- c(
  list.files(input_dir, full.names = TRUE),
  list.files(expected_dir, full.names = TRUE),
  file.path(fixture_root, "analysis-parameters.json"),
  file.path(fixture_root, "known-answer-claims.csv"),
  file.path(fixture_root, "README.md")
)
inventory <- data.frame(
  relative_path = substring(inventory_paths, nchar(fixture_root) + 2L),
  size_bytes = as.numeric(file.info(inventory_paths)$size),
  md5 = unname(tools::md5sum(inventory_paths)),
  stringsAsFactors = FALSE
)
inventory <- inventory[order(inventory$relative_path), , drop = FALSE]
utils::write.csv(
  inventory,
  file.path(fixture_root, "checksums.csv"),
  row.names = FALSE
)
canonicalize_lf(file.path(fixture_root, "checksums.csv"))

message(
  "Generated ", fixture_root, " with ",
  length(tracks), " streamlines and ", nrow(inventory),
  " checksummed files."
)
