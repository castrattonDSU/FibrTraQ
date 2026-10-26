test_that("output guide preserves existing guides and follows its own portable links", {
  root <- tempfile("output-guide- # 100%20 ")
  dir.create(root)
  first <- FibrTraQ:::ftq_write_output_guide(root)
  expect_identical(basename(first), "OUTPUTS_README.html")
  before <- tools::md5sum(list.files(root, full.names = TRUE))
  expect_identical(FibrTraQ:::ftq_write_output_guide(root), first)
  expect_identical(tools::md5sum(names(before)), before)
  expect_match(paste(readLines(first, encoding = "UTF-8"), collapse = "\n"),
               'href="OUTPUTS_README.md"', fixed = TRUE)
  expect_true(file.exists(file.path(dirname(first), "OUTPUTS_README.md")))

  # A preserved old/custom guide must not prevent a newer report from rendering.
  writeLines("Original investigator notes", first)
  old_hash <- tools::md5sum(first)
  versioned <- FibrTraQ:::ftq_write_output_guide(root)
  expect_match(basename(versioned), "^OUTPUTS_README-[a-f0-9]{12}[.]html$")
  expect_identical(tools::md5sum(first), old_hash)
  expect_identical(FibrTraQ:::ftq_write_output_guide(root), versioned)
  md <- sub("[.]html$", ".md", basename(versioned))
  expect_true(file.exists(file.path(root, md)))
  expect_match(paste(readLines(versioned, encoding = "UTF-8"), collapse = "\n"),
               paste0('href="', md, '"'), fixed = TRUE)
  expect_equal(length(list.files(root)), 4L)

  writeLines("Further investigator notes", versioned)
  all_hashes <- tools::md5sum(list.files(root, full.names = TRUE))
  expect_error(FibrTraQ:::ftq_write_output_guide(root), "Content-specific output guide")
  expect_identical(tools::md5sum(names(all_hashes)), all_hashes)
})

test_that("dictionary roles preserve target/coordinate collisions and reject undocumented fields", {
  definitions <- FibrTraQ:::ftq_output_field_dictionary()
  expect_false(anyDuplicated(paste(definitions$family, definitions$column)) > 0L)
  expect_true(all(vapply(definitions, function(x) all(!is.na(x) & nzchar(x)), logical(1))))
  expect_setequal(unique(definitions$family), names(FibrTraQ:::ftq_output_family_files()))
  expect_error(FibrTraQ:::ftq_output_field_definition("unknown", "x"), "Unknown output family")
  expect_error(FibrTraQ:::ftq_output_field_definition("per_voxel", "unexpected"), "Undocumented")
  source <- data.frame(fNode_Coords = c("1_1_1", "2_1_1"), tracks = c(2L, 0L),
                       local = c(1L, 0L), leave_no_target = c(0L, 0L),
                       x = c(1L, 0L), voxel_x = c(1L, 0L), none = c(0L, 1L))
  view <- FibrTraQ:::ftq_report_voxel_export(source)
  coordinates <- view$columns$column[view$columns$role == "coordinate"]
  targets <- view$columns$column[view$columns$role == "target"]
  covered <- FibrTraQ:::ftq_output_field_definition("all_voxels", names(view$data), targets, coordinates)
  expect_identical(covered$column, names(view$data))
  expect_true(all(covered$units[covered$column %in% coordinates] == "one-based voxel index"))
  expect_true(all(covered$units[covered$column %in% targets] == "count"))
  expect_error(FibrTraQ:::ftq_output_field_definition("per_voxel", names(view$data), targets), "Undocumented")
  expect_silent(FibrTraQ:::ftq_output_field_definition("all_voxels_columns", names(view$columns)))
  expect_silent(FibrTraQ:::ftq_output_field_definition("local_fraction", names(FibrTraQ:::ftq_report_local_fraction(source))))
  preview <- FibrTraQ:::ftq_output_field_definition("top_voxels",
    c("voxel_x_2", "y", "z", "leave_with_target", "derived_leave_with_target"),
    target_labels = "leave_with_target", coordinate_columns = c("voxel_x_2", "y", "z"),
    derived_column = "derived_leave_with_target")
  expect_match(preview$definition[preview$column == "leave_with_target"], "named target", fixed = TRUE)
  expect_match(preview$definition[preview$column == "derived_leave_with_target"], "tracks - local - leave_no_target", fixed = TRUE)
  aliases <- definitions[definitions$family == "rsc", ]
  expect_equal(aliases$alias_of[aliases$column == "tracks"], "target_associations")
  expect_equal(aliases$alias_of[aliases$column == "total_origin_voxels"], "total_focal_voxels")
})

test_that("dictionary covers generated scientific and supporting table schemas", {
  blank <- array(0, dim = c(8L, 4L, 4L))
  focal <- target <- blank
  focal[2:3, 2, 2] <- 1
  target[6, 2, 2] <- 1
  attr(focal, "fibrtraq_voxel_to_world_basis") <- diag(3)
  paths <- list(rbind(c(2, 2, 2), c(3, 2, 2), c(6, 2, 2)),
                rbind(c(1, 2, 2), c(2, 2, 2), c(6, 2, 2), c(7, 2, 2)),
                matrix(c(2, 2, 2), nrow = 1L))
  check <- function(family, table, targets = "T", coords = character()) {
    expect_identical(FibrTraQ:::ftq_output_field_definition(family, names(table), targets, coords)$column,
                     names(table), info = family)
  }
  for (mode in c("full_streamline", "endpoint_voxels", "local_search")) {
    result <- FibrTraQ::FibrTraQ(fNode = focal, cNode_list = list(T = target), cNode_labels = "T",
      localROI = array(0, dim = dim(focal)), all_tracks = paths, track_space = "voxel",
      discretization = "round", axis_flip = c(FALSE, FALSE, FALSE), backend = "r",
      assignment_mode = mode, search_radius_mm = 1, endpoint_analysis = "exact_and_dilated1",
      endpoint_context_list = list(brain = array(1, dim = dim(focal))))
    roi <- attr(result, "fibrtraq_roi_summary")
    rsc <- FibrTraQ::rsc_summary(result, percent = TRUE, warn_legacy = FALSE)
    endpoints <- FibrTraQ::ftq_endpoint_tables(result)
    check("per_voxel", result)
    check("roi_summary", roi)
    check("rsc", rsc)
    check("fate", FibrTraQ:::ftq_fate_summary(result))
    check("connection_summary", FibrTraQ:::ftq_report_connection_summary(result, rsc, roi))
    check("report_fate", FibrTraQ:::ftq_report_fate_summary(result, roi))
    check("endpoint_summary", endpoints$summary)
    check("endpoint_penetration", endpoints$penetration)
    check("endpoint_context", endpoints$context)
    check("penetration_summary", FibrTraQ:::ftq_endpoint_penetration_summary(endpoints$penetration))
    top_dir <- tempfile("guide-preview-")
    dir.create(top_dir)
    tops <- FibrTraQ:::ftq_write_top_tables(result, top_dir, "test", n = 50L)
    for (path in tops) check("top_voxels", utils::read.csv(path, check.names = FALSE))
    exploratory <- FibrTraQ::ftq_exploratory_summary(result, write_maps = FALSE)
    check("exploratory_per_voxel", exploratory$data)
    check("exploratory_target_totals", exploratory$target_totals)
    check("exploratory_lut", FibrTraQ:::ftq_exploratory_lut("T"))
    check("exploratory_provenance", exploratory$provenance)
    summary <- FibrTraQ::ftq_stat_subject_summary(result, endpoint_summary = endpoints$summary,
                                                endpoint_penetration = endpoints$penetration)
    for (nm in names(summary)) check(paste0("stat_", nm), summary[[nm]])
    target_tables <- FibrTraQ::ftq_target_summary_tables(result, roi, endpoint_summary = endpoints$summary)
    check("target_unique", target_tables$unique_streamlines)
    check("target_associations", target_tables$voxel_associations)
    check("target_endpoints", target_tables$endpoint_streamlines)
    check("spatial_gradients", FibrTraQ::ftq_stat_spatial_gradients(result))
    row <- FibrTraQ:::ftq_regional_study_run_row("001", "focal", NULL, NULL, NULL, "test", NULL)
    check("run_index", row)
    check("validation_suite", FibrTraQ:::ftq_validation_suite_row(manifest_path = "test", run_status = "not_run", started_at = "test"))
    if (mode == "local_search") {
      local <- FibrTraQ::ftq_local_search_tables(result)
      check("local_search_assignments", local$assignments)
      check("local_search_summary", local$summary)
      check("local_search_endpoint_status", local$endpoint_status)
    }
  }
})

test_that("installed output guide is synchronized with its executable definitions", {
  asset <- system.file("output-guide", "OUTPUTS_README.md", package = "FibrTraQ")
  expect_true(nzchar(asset))
  expect_identical(readLines(asset, encoding = "UTF-8"), FibrTraQ:::ftq_output_guide_markdown())
})
