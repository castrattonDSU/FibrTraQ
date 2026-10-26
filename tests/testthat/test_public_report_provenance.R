test_that("primary report manifest recovers exact companion identity and method", {
  root <- tempfile("primary-provenance-")
  dir.create(root)
  input <- file.path(root, "001_fibrtraq_per-voxel-hits.csv")
  metadata <- file.path(root, "001_fibrtraq_run-metadata.json")
  writeLines("Not read by the metadata resolver", input)
  saved <- list(subject = "001", prefix = "001-cerebellum", fnode = "C:/original/focal.nii.gz",
    target_intersection_mode = "local_search", local_search_radius_mm = 3.25,
    focal_assignment_mode = "first_intersection", endpoint_analysis = "exact_and_dilated1")
  jsonlite::write_json(list(manifest = saved,
    summary = list(target_intersection_mode = "local_search", local_search_radius_mm = 3.25),
    files = list(result_csv = "C:\\original\\001_fibrtraq_per-voxel-hits.csv")),
    metadata, auto_unbox = TRUE)
  resolved <- FibrTraQ:::ftq_report_resolve_manifest(input, NULL)
  expect_identical(resolved$subject, "001")
  expect_identical(resolved$prefix, "001-cerebellum")
  expect_identical(resolved$fnode, "C:/original/focal.nii.gz")
  expect_identical(resolved$target_intersection_mode, "local_search")
  expect_identical(resolved$local_search_radius_mm, 3.25)
  expect_identical(resolved$focal_assignment_mode, "first_intersection")
  expect_null(resolved$provenance_warning)
  expect_identical(resolved$provenance_metadata_json, metadata)
  explicit <- list(subject = "S1c", prefix = "new-view", analysis = list(
    target_intersection_mode = "local_search", local_search_radius_mm = 3.25))
  combined <- FibrTraQ:::ftq_report_resolve_manifest(input, explicit)
  expect_identical(combined$subject, "S1c")
  expect_identical(combined$prefix, "new-view")
  expect_identical(combined$analysis$target_intersection_mode, "local_search")
  expect_error(FibrTraQ:::ftq_report_resolve_manifest(input,
    list(target_intersection_mode = "endpoint_voxels")), "conflicts with recorded")
  expect_error(FibrTraQ:::ftq_report_resolve_manifest(input,
    list(target_intersection_mode = "local_search", local_search_radius_mm = 2)),
    "conflicts with recorded")
  saved$analysis <- list(target_intersection_mode = "local_search", endpoint_analysis = "exact_and_dilated1")
  jsonlite::write_json(list(manifest = saved,
    files = list(result_csv = input)), metadata, auto_unbox = TRUE)
  endpoint_override <- FibrTraQ:::ftq_report_resolve_manifest(input,
    list(target_intersection_mode = "local_search", endpoint_analysis = "none"))
  expect_identical(endpoint_override$endpoint_analysis, "none")
  expect_null(endpoint_override$analysis$endpoint_analysis)
  jsonlite::write_json(list(manifest = saved,
    summary = list(target_intersection_mode = "endpoint_voxels")), metadata, auto_unbox = TRUE)
  expect_error(FibrTraQ:::ftq_report_resolve_manifest(input, NULL), "Conflicting target_intersection_mode")
  jsonlite::write_json(list(files = list(result_csv = "other.csv")), metadata, auto_unbox = TRUE)
  expect_error(FibrTraQ:::ftq_report_resolve_manifest(input, NULL), "different result CSV")
})

test_that("primary report uses result attributes and canonical assignment aliases", {
  result <- data.frame()
  for (mode in c("full_streamline", "endpoint_voxels", "local_search")) {
    attr(result, "fibrtraq_target_intersection_mode") <- mode
    attr(result, "fibrtraq_local_search_radius_mm") <- 4.5
    attr(result, "fibrtraq_focal_assignment_mode") <- "first_intersection"
    resolved <- FibrTraQ:::ftq_report_resolve_manifest(result, NULL)
    expect_identical(resolved$target_intersection_mode, mode)
    expect_identical(resolved$focal_assignment_mode, "first_intersection")
    expect_null(resolved$provenance_warning)
    if (mode == "local_search") expect_identical(resolved$local_search_radius_mm, 4.5)
  }
  attr(result, "fibrtraq_target_intersection_mode") <- "endpoint_voxels"
  resolved <- FibrTraQ:::ftq_report_resolve_manifest(result,
    list(analysis = list(target_intersection_mode = "end_voxels")))
  expect_identical(resolved$target_intersection_mode, "endpoint_voxels")
  expect_identical(resolved$analysis$target_intersection_mode, "endpoint_voxels")
})

test_that("unknown primary method and radius are visibly distinguished from endpoint QC", {
  expect_warning(resolved <- FibrTraQ:::ftq_report_resolve_manifest("unknown.csv",
    list(endpoint_analysis = "exact_and_dilated1")), "not a verified method")
  expect_identical(resolved$target_intersection_mode, "full_streamline")
  expect_match(resolved$provenance_warning, "Endpoint QC settings do not establish", fixed = TRUE)
  expect_identical(resolved$endpoint_analysis, "exact_and_dilated1")
  expect_warning(local <- FibrTraQ:::ftq_report_resolve_manifest("unknown.csv",
    list(analysis = list(target_intersection_mode = "local_search"))), "2 mm report default is unverified")
  expect_identical(local$local_search_radius_mm, 2)
  expect_identical(local$analysis$local_search_radius_mm, 2)
  expect_match(local$provenance_warning, "Local-search radius was not recorded", fixed = TRUE)
  # Legacy untagged data frames remain quiet but the HTML receives the warning.
  expect_warning(frame <- FibrTraQ:::ftq_report_resolve_manifest(data.frame(), NULL), NA)
  expect_match(frame$provenance_warning, "not a verified method", fixed = TRUE)
  expect_warning(known <- FibrTraQ:::ftq_report_resolve_manifest("unknown.csv",
    list(target_intersection_mode = "full_streamline")), NA)
  expect_null(known$provenance_warning)
})

test_that("manifest-file defaults do not erase missing primary method provenance", {
  root <- tempfile("manifest-file-provenance-")
  dir.create(root)
  manifest_path <- file.path(root, "manifest.json")
  manifest <- list(subject = "001", description = "original", dataset_root = root,
    fnode = "focal.nii.gz", tractogram = "streamlines.trk",
    targets = list(list(label = "A", path = "target.nii.gz")),
    analysis = list(endpoint_analysis = "none"))
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE)
  expect_warning(resolved <- FibrTraQ:::ftq_report_resolve_manifest("unknown.csv", manifest_path),
    "not a verified method")
  expect_identical(resolved$subject, "001")
  expect_identical(
    chartr("\\", "/", resolved$fnode),
    file.path(normalizePath(root, winslash = "/"), "focal.nii.gz")
  )
  expect_identical(resolved$analysis$endpoint_analysis, "none")
  expect_identical(resolved$target_intersection_mode, "full_streamline")
  expect_match(resolved$provenance_warning, "not a verified method", fixed = TRUE)
})

test_that("legacy writer NA metadata is missing evidence rather than a method or radius", {
  root <- tempfile("legacy-missing-provenance-")
  result <- data.frame(fNode_Coords = c("1_1_1", "2_2_2"),
    tracks = c(2L, 0L), local = c(1L, 0L), leave_no_target = c(0L, 0L),
    A = c(1L, 0L), none = c(0L, 1L))
  out <- ftq_write_outputs(result, rsc_summary(result, warn_legacy = FALSE),
    template_image = NULL, output_dir = root, prefix = "untagged", write_maps = FALSE)
  # Missing legacy attributes can be JSON null or the text "NA", depending on
  # their R type. Exercise the explicit sentinel on an actual writer sidecar.
  metadata <- jsonlite::read_json(out$metadata_json)
  metadata$summary$target_intersection_mode <- "NA"
  metadata$summary$local_search_radius_mm <- "NA"
  jsonlite::write_json(metadata, out$metadata_json, auto_unbox = TRUE)
  expect_warning(resolved <- FibrTraQ:::ftq_report_resolve_manifest(out$result_csv, NULL),
    "not a verified method")
  expect_identical(resolved$target_intersection_mode, "full_streamline")
  expect_warning(local <- FibrTraQ:::ftq_report_resolve_manifest(out$result_csv,
    list(analysis = list(target_intersection_mode = "local_search", local_search_radius_mm = 9))), NA)
  expect_identical(local$local_search_radius_mm, 9)
  expect_warning(exploratory <- ftq_exploratory_summary(out$result_csv, write_maps = FALSE),
    "not a verified method")
  expect_identical(exploratory$provenance$target_intersection_mode, "full_streamline")
  metadata$summary$target_intersection_mode <- "unrecognized_method"
  jsonlite::write_json(metadata, out$metadata_json, auto_unbox = TRUE)
  expect_error(FibrTraQ:::ftq_report_resolve_manifest(out$result_csv, NULL), "must be one of")
})

test_that("primary HTML exposes unverified method and rejects contradictory CSV relabeling", {
  root <- tempfile("primary-method-html-")
  result <- data.frame(fNode_Coords = c("1_1_1", "2_2_2"),
    tracks = c(2L, 0L), local = c(1L, 0L), leave_no_target = c(0L, 0L),
    A = c(1L, 0L), none = c(0L, 1L))
  rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
  expect_warning(report <- ftq_report(result, rsc,
    output_dir = file.path(root, "unverified")), NA)
  html <- paste(readLines(report$html), collapse = "")
  expect_match(html, "not a verified method for these data", fixed = TRUE)
  expect_match(html, "Endpoint QC settings do not establish", fixed = TRUE)
  attr(result, "fibrtraq_target_intersection_mode") <- "full_streamline"
  out <- ftq_write_outputs(result, rsc, template_image = NULL,
    output_dir = file.path(root, "saved"), prefix = "known", write_maps = FALSE)
  expect_error(ftq_report(out$result_csv, out$rsc_csv,
    manifest = list(target_intersection_mode = "local_search", local_search_radius_mm = 3.25),
    output_dir = file.path(root, "contradiction")), "conflicts with recorded")
  expect_false(dir.exists(file.path(root, "contradiction")))
})

test_that("focal allocation and voxelization cannot contradict recorded report methods", {
  choices <- list(
    focal_assignment_mode = c("first_intersection", "all_focal_voxels"),
    streamline_voxelization_mode = c("segment_interpolation", "sample_points")
  )
  root <- tempfile("additional-method-provenance-")
  dir.create(root)
  input <- file.path(root, "known_fibrtraq_per-voxel-hits.csv")
  metadata_path <- file.path(root, "known_fibrtraq_run-metadata.json")
  for (field in names(choices)) {
    method <- choices[[field]][[1L]]
    contradictory <- choices[[field]][[2L]]
    result <- data.frame()
    attr(result, "fibrtraq_target_intersection_mode") <- "full_streamline"
    attr(result, paste0("fibrtraq_", field)) <- method
    explicit <- list(target_intersection_mode = "full_streamline")
    explicit[[field]] <- contradictory
    expect_error(FibrTraQ:::ftq_report_resolve_manifest(result, explicit),
      paste0(field, " conflicts with recorded"))
    # Partial names accepted by the existing manifest interface are normalized
    # before comparison and written back as canonical method values.
    explicit[[field]] <- substr(method, 1L, 5L)
    normalized <- FibrTraQ:::ftq_report_resolve_manifest(result, explicit)
    expect_identical(normalized[[field]], method)

    manifest <- list(target_intersection_mode = "full_streamline")
    manifest[[field]] <- method
    saved_summary <- list(target_intersection_mode = "full_streamline")
    saved_summary[[field]] <- method
    metadata <- list(manifest = manifest, summary = saved_summary)
    jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE)
    explicit <- list(analysis = stats::setNames(list(contradictory), field))
    expect_error(FibrTraQ:::ftq_report_resolve_manifest(input, explicit),
      paste0(field, " conflicts with recorded"))
    metadata$manifest$analysis <- stats::setNames(list(contradictory), field)
    jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE)
    expect_error(FibrTraQ:::ftq_report_resolve_manifest(input, NULL),
      paste0("Conflicting ", field))
    metadata$manifest$analysis <- NULL
    metadata$summary[[field]] <- contradictory
    jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE)
    expect_error(FibrTraQ:::ftq_report_resolve_manifest(input, NULL),
      paste0("Conflicting ", field))
    metadata$summary[[field]] <- "invalid_method"
    jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE)
    expect_error(FibrTraQ:::ftq_report_resolve_manifest(input, NULL), "should be one of")
    metadata$summary[[field]] <- "NA"
    metadata$manifest[[field]] <- "NA"
    jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE)
    expect_warning(legacy <- FibrTraQ:::ftq_report_resolve_manifest(input, NULL), NA)
    expect_null(legacy[[field]])
  }
})
