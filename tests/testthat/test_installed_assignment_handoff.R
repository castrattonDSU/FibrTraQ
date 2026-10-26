test_that("installed file workflow records endpoint and local assignment modes", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("freesurferformats")

  public_arguments <- names(formals(ftq_run_files))
  expect_true("assignment_mode" %in% public_arguments)
  expect_true("search_radius_mm" %in% public_arguments)

  example <- ftq_example_data(verify = TRUE)
  params <- example$parameters
  modes <- data.frame(
    mode = c("endpoint_voxels", "local_search"),
    radius_mm = c(2, 1.5),
    stringsAsFactors = FALSE
  )
  run_mode <- function(mode, radius_mm, output_dir) {
    ftq_run_files(
      tractogram = example$tractogram,
      focal_mask = example$focal_mask,
      target_masks = example$target_masks,
      output_dir = output_dir,
      subject = "handoff",
      focal_label = "focal",
      local_mask = example$local_mask,
      reference_image = example$reference_image,
      description = paste0("handoff-", mode),
      track_space = params$track_space,
      discretization = params$discretization,
      axis_flip = as.logical(params$axis_flip),
      focal_assignment_mode = params$focal_assignment_mode,
      streamline_voxelization_mode = params$streamline_voxelization_mode,
      endpoint_analysis = "none",
      track_reader = "memory",
      backend = "r",
      write_maps = FALSE,
      write_report = TRUE,
      write_exploratory_report = FALSE,
      write_qc_viewer = FALSE,
      overwrite = TRUE,
      verbose = FALSE,
      progress = FALSE,
      assignment_mode = mode,
      search_radius_mm = radius_mm
    )
  }

  for (i in seq_len(nrow(modes))) {
    mode <- modes$mode[[i]]
    radius_mm <- modes$radius_mm[[i]]
    output_dir <- tempfile(paste0("fibrtraq-handoff-", mode, "-"))
    on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

    if (identical(mode, "local_search")) {
      expect_warning(
        result <- run_mode(mode, radius_mm, output_dir),
        "does not declare NIfTI spatial units",
        fixed = TRUE
      )
    } else {
      expect_no_warning(result <- run_mode(mode, radius_mm, output_dir))
    }

    expect_s3_class(result, "ftq_regional_study_result")
    expect_equal(nrow(result$index), 1L)
    run <- result$index[1L, , drop = FALSE]
    expect_identical(run$run_status, "completed")
    expect_true(run$validation_passed)
    expect_identical(run$target_intersection_mode, mode)
    expect_equal(run$local_search_radius_mm, radius_mm)
    expect_identical(run$error, "")

    artifacts <- unname(c(
      result$index_csv,
      result$report_index,
      run$result_csv,
      run$rsc_csv,
      run$fate_csv,
      run$roi_summary_csv,
      run$metadata_json,
      run$validation_json,
      run$primary_report
    ))
    expect_true(all(file.exists(artifacts)))
    expect_true(all(file.info(artifacts)$size > 0))

    run_index <- utils::read.csv(
      result$index_csv,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    roi <- utils::read.csv(
      run$roi_summary_csv,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    voxel <- utils::read.csv(
      run$result_csv,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    validation <- jsonlite::read_json(
      run$validation_json,
      simplifyVector = TRUE
    )
    metadata <- jsonlite::read_json(
      run$metadata_json,
      simplifyVector = TRUE
    )

    expect_identical(run_index$target_intersection_mode, mode)
    expect_equal(run_index$local_search_radius_mm, radius_mm)
    expect_identical(roi$target_intersection_mode, mode)
    expect_true(validation$passed)
    expect_identical(
      metadata$manifest$analysis$target_intersection_mode,
      mode
    )
    expect_identical(metadata$manifest$target_intersection_mode, mode)
    expect_identical(metadata$summary$target_intersection_mode, mode)
    expect_equal(metadata$manifest$analysis$local_search_radius_mm, radius_mm)
    expect_equal(metadata$manifest$local_search_radius_mm, radius_mm)
    if (identical(mode, "local_search")) {
      expect_equal(roi$local_search_radius_mm, radius_mm)
      expect_equal(metadata$summary$local_search_radius_mm, radius_mm)
    } else {
      inactive_roi_radius <- roi[["local_search_radius_mm"]]
      inactive_summary_radius <- suppressWarnings(as.numeric(
        metadata$summary$local_search_radius_mm
      ))
      expect_true(!length(inactive_roi_radius) || all(is.na(inactive_roi_radius)))
      expect_true(!length(inactive_summary_radius) ||
                    all(is.na(inactive_summary_radius)))
    }

    expect_true(all(voxel$tracks >= 0))
    expect_true(all(voxel$local >= 0))
    expect_true(all(voxel$leave_no_target >= 0))
    expect_true(all(voxel$local + voxel$leave_no_target <= voxel$tracks))
    target_columns <- intersect(names(example$target_masks), names(voxel))
    expect_setequal(target_columns, names(example$target_masks))
    expect_true(all(vapply(target_columns, function(target) {
      all(voxel[[target]] >= 0 & voxel[[target]] <= voxel$tracks)
    }, logical(1))))
    expect_lte(
      roi$unique_any_target_positive_streamlines,
      roi$unique_focal_intersecting_streamlines
    )
    expect_true(roi$unique_exclusive_fate_sum_matches_focal)

    report_html <- paste(
      readLines(run$primary_report, warn = FALSE),
      collapse = "\n"
    )
    expect_match(
      report_html,
      paste0("<b>Target assignment mode:</b> ", mode),
      fixed = TRUE
    )
    if (identical(mode, "local_search")) {
      expect_match(
        report_html,
        paste0("<b>Local search radius (mm):</b> ", radius_mm),
        fixed = TRUE
      )
    } else {
      expect_false(grepl("<b>Local search radius (mm):</b>", report_html,
                         fixed = TRUE))
    }
  }
})
