library(testthat)
library(FibrTraQ)

ftq_target_assignment_mask <- function(dims, x) {
  out <- array(0, dim = dims)
  if (length(x)) out[x, 1, 1] <- 1
  out
}

ftq_target_assignment_point_mask <- function(dims, points) {
  out <- array(0, dim = dims)
  points <- as.matrix(points)
  if (nrow(points)) out[points] <- 1
  out
}

ftq_target_assignment_run <- function(focal,
                                      targets,
                                      tracks,
                                      mode = NULL,
                                      backend = "r",
                                      radius_mm = 2,
                                      voxel_to_world = diag(3),
                                      discretization = "round",
                                      voxelization = "sample_points",
                                      reverse_streamlines = FALSE,
                                      endpoint_analysis = "none") {
  attr(focal, "fibrtraq_voxel_to_world_basis") <-
    as.matrix(voxel_to_world)
  args <- list(
    fNode = focal,
    localROI = array(0, dim = dim(focal)),
    cNode_list = targets,
    cNode_labels = names(targets),
    all_tracks = tracks,
    discretization = discretization,
    focal_assignment_mode = "all_focal_voxels",
    streamline_voxelization_mode = voxelization,
    reverse_streamlines = reverse_streamlines,
    track_space = "voxel",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = backend,
    endpoint_analysis = endpoint_analysis,
    search_radius_mm = radius_mm
  )
  if (!is.null(mode)) args$assignment_mode <- mode
  do.call(FibrTraQ, args)
}

ftq_target_assignment_core <- function(x) {
  attributes(x) <- attributes(x)[c("names", "row.names", "class")]
  x
}

test_that("target assignment modes are selectable and aliases preserve the default", {
  dims <- c(10L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 6L)
  tracks <- list(
    rbind(c(2, 1, 1), c(4, 1, 1), c(6, 1, 1)),
    rbind(c(1, 1, 1), c(2, 1, 1), c(6, 1, 1), c(7, 1, 1)),
    rbind(c(1, 1, 1), c(2, 1, 1), c(6, 1, 1), c(9, 1, 1))
  )

  backend_results <- lapply(c("r", "cpp"), function(backend) {
    default <- ftq_target_assignment_run(
      focal, list(T = target), tracks, backend = backend
    )
    full <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "full_streamline", backend = backend
    )
    all_alias <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "all_voxels", backend = backend
    )
    endpoint <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "endpoint_voxels", backend = backend
    )
    endpoint_alias <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "end_voxels", backend = backend
    )
    local <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "local_search", backend = backend
    )

    expect_equal(default$tracks, 3L, info = backend)
    expect_equal(c(full$T, endpoint$T, local$T), c(3L, 1L, 2L),
                 info = backend)
    expect_equal(
      ftq_target_assignment_core(default),
      ftq_target_assignment_core(full),
      info = backend
    )
    expect_equal(
      ftq_target_assignment_core(all_alias),
      ftq_target_assignment_core(full),
      info = backend
    )
    expect_equal(
      ftq_target_assignment_core(endpoint_alias),
      ftq_target_assignment_core(endpoint),
      info = backend
    )
    expect_equal(
      attr(all_alias, "fibrtraq_target_intersection_mode"),
      "full_streamline",
      info = backend
    )
    expect_equal(
      attr(endpoint_alias, "fibrtraq_target_intersection_mode"),
      "endpoint_voxels",
      info = backend
    )
    expect_equal(
      attr(local, "fibrtraq_local_search_radius_mm"),
      2,
      info = backend
    )
    list(full = full, endpoint = endpoint, local = local)
  })

  for (mode in c("full", "endpoint", "local")) {
    expect_equal(
      ftq_target_assignment_core(backend_results[[1L]][[mode]]),
      ftq_target_assignment_core(backend_results[[2L]][[mode]]),
      info = mode
    )
  }
})

test_that("local-search audit reconciles proximity-only and full-path connections", {
  dims <- c(12L, 2L, 2L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 10L)
  tracks <- list(
    rbind(c(2, 1, 1), c(10, 1, 1)),
    rbind(c(1.4, 1, 1), c(9.4, 1, 1)),
    rbind(c(2, 1, 1), c(7, 1, 1))
  )

  results <- lapply(c("r", "cpp"), function(backend) {
    result <- ftq_target_assignment_run(
      focal,
      list(T = target),
      tracks,
      mode = "local_search",
      backend = backend,
      radius_mm = 0.75,
      voxelization = "segment_interpolation",
      endpoint_analysis = "none"
    )
    audit <- ftq_local_search_tables(result)
    roi <- attr(result, "fibrtraq_roi_summary")
    checks <- FibrTraQ:::ftq_validate_local_search_assignments(
      audit$assignments,
      roi_summary = roi,
      target_labels = "T"
    )

    expect_true(audit$metadata$available, info = backend)
    expect_equal(nrow(audit$assignments), 3L, info = backend)
    expect_equal(result$T, 2L, info = backend)
    expect_equal(
      audit$summary$connections_with_full_path_intersection,
      1L,
      info = backend
    )
    expect_equal(
      audit$summary$connections_without_full_path_intersection,
      1L,
      info = backend
    )
    expect_equal(
      audit$summary$full_path_intersection_percent_of_connections,
      50,
      info = backend
    )
    expect_equal(
      audit$summary$full_path_target_distance_median_mm,
      0,
      tolerance = 1e-10,
      info = backend
    )
    expect_equal(
      audit$summary$proximity_only_target_distance_median_mm,
      0.6,
      tolerance = 1e-10,
      info = backend
    )
    expect_equal(
      sort(audit$assignments$paired_max_distance_mm[!is.na(
        audit$assignments$paired_target
      )]),
      c(0, 0.6),
      tolerance = 1e-10,
      info = backend
    )
    expect_false(any(grepl(
      "unique_target_streamlines_(with|without)_full_path_intersection",
      names(roi)
    )), info = backend)
    expect_true(all(checks$passed), info = backend)
    list(result = result, audit = audit, roi = roi)
  })

  expect_equal(
    ftq_target_assignment_core(results[[1L]]$result),
    ftq_target_assignment_core(results[[2L]]$result)
  )
  expect_equal(
    results[[1L]]$audit,
    results[[2L]]$audit,
    tolerance = 1e-10
  )

  tampered <- results[[1L]]$audit$assignments
  paired_row <- which(!is.na(tampered$paired_target))[[1L]]
  tampered$paired_target_intersects_full_path[paired_row] <- NA
  tampered_checks <- FibrTraQ:::ftq_validate_local_search_assignments(
    tampered,
    roi_summary = results[[1L]]$roi,
    target_labels = "T"
  )
  expect_false(tampered_checks$passed[
    tampered_checks$check ==
      "local_search_full_path_partition_reconciles:T"
  ])
})

test_that("local-search audit is invariant to streamline reversal", {
  dims <- c(12L, 2L, 2L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 10L)
  tracks <- list(
    rbind(c(2, 1, 1), c(10, 1, 1)),
    rbind(c(1.4, 1, 1), c(9.4, 1, 1))
  )

  for (backend in c("r", "cpp")) {
    forward <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "local_search", backend = backend,
      radius_mm = 0.75, voxelization = "segment_interpolation"
    )
    reversed <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "local_search", backend = backend,
      radius_mm = 0.75, voxelization = "segment_interpolation",
      reverse_streamlines = TRUE
    )
    forward_audit <- ftq_local_search_tables(forward)
    reversed_audit <- ftq_local_search_tables(reversed)

    expect_equal(
      ftq_target_assignment_core(forward),
      ftq_target_assignment_core(reversed),
      info = backend
    )
    expect_equal(
      forward_audit$summary,
      reversed_audit$summary,
      tolerance = 1e-10,
      info = backend
    )
    expect_equal(
      forward_audit$assignments$start_status,
      reversed_audit$assignments$end_status,
      info = backend
    )
    expect_equal(
      forward_audit$assignments$end_status,
      reversed_audit$assignments$start_status,
      info = backend
    )
    expect_equal(
      forward_audit$assignments$start_distance_mm,
      reversed_audit$assignments$end_distance_mm,
      tolerance = 1e-10,
      info = backend
    )
    expect_equal(
      forward_audit$assignments$end_distance_mm,
      reversed_audit$assignments$start_distance_mm,
      tolerance = 1e-10,
      info = backend
    )
  }
})

test_that("non-local results omit the additive local-search audit", {
  dims <- c(8L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 6L)
  result <- ftq_target_assignment_run(
    focal,
    list(T = target),
    list(rbind(c(2, 1, 1), c(6, 1, 1))),
    mode = "full_streamline"
  )
  audit <- ftq_local_search_tables(result)

  expect_false(audit$metadata$available)
  expect_equal(nrow(audit$assignments), 0L)
  expect_equal(nrow(audit$summary), 0L)
})

test_that("endpoint-based modes preserve hand-calculated focal and fate counts", {
  dims <- c(10L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, c(3L, 4L))
  local <- ftq_target_assignment_mask(dims, c(3L, 4L, 6L, 8L))
  target <- ftq_target_assignment_mask(dims, 8L)
  attr(focal, "fibrtraq_voxel_to_world_basis") <- diag(3)

  tracks <- list(
    # Local/positive; revisits x=3 and occupies both focal voxels.
    rbind(c(3, 1, 1), c(4, 1, 1), c(3, 1, 1), c(8, 1, 1)),
    # Local/negative; revisits one focal voxel.
    rbind(c(3, 1, 1), c(3, 1, 1), c(6, 1, 1)),
    # Nonlocal/positive; revisits the other focal voxel.
    rbind(c(4, 1, 1), c(4, 1, 1), c(1, 1, 1), c(8, 1, 1)),
    # Nonlocal/negative; occupies both focal voxels.
    rbind(c(4, 1, 1), c(3, 1, 1), c(1, 1, 1), c(6, 1, 1))
  )
  expected_result <- data.frame(
    fNode_Coords = c("3_1_1", "4_1_1"),
    tracks = c(3L, 3L),
    local = c(2L, 1L),
    leave_no_target = c(1L, 1L),
    T = c(1L, 2L),
    none = c(0L, 0L),
    stringsAsFactors = FALSE
  )
  expected_roi <- c(
    total_input_streamlines = 4,
    processed_streamlines = 4,
    unique_focal_intersecting_streamlines = 4,
    unique_local_streamlines = 2,
    unique_leave_no_target_streamlines = 1,
    unique_leave_with_target_streamlines = 1,
    unique_focal_intersecting_with_target_streamlines = 2,
    unique_local_target_positive_streamlines = 1,
    unique_local_target_negative_streamlines = 1,
    unique_nonlocal_target_positive_streamlines = 1,
    unique_nonlocal_target_negative_streamlines = 1,
    unique_any_target_positive_streamlines = 2,
    unique_any_target_negative_streamlines = 2,
    unique_exclusive_fate_sum = 4,
    total_focal_voxel_streamline_associations = 6,
    total_focal_voxel_streamline_target_associations = 3,
    `unique_target_streamlines:T` = 2,
    `unique_full_path_target_streamlines:T` = 2
  )

  backend_results <- list()
  for (mode in c("endpoint_voxels", "local_search")) {
    backend_results[[mode]] <- list()
    for (backend in c("r", "cpp")) {
      result <- FibrTraQ(
        fNode = focal,
        localROI = local,
        cNode_list = list(T = target),
        cNode_labels = "T",
        all_tracks = tracks,
        discretization = "round",
        focal_assignment_mode = "all_focal_voxels",
        streamline_voxelization_mode = "sample_points",
        reverse_streamlines = FALSE,
        track_space = "voxel",
        axis_flip = c(FALSE, FALSE, FALSE),
        backend = backend,
        endpoint_analysis = "none",
        assignment_mode = mode,
        search_radius_mm = 0.25
      )
      roi <- attr(result, "fibrtraq_roi_summary")

      expect_equal(
        result[names(expected_result)],
        expected_result,
        info = paste(mode, backend, "per-voxel")
      )
      for (field in names(expected_roi)) {
        expect_equal(
          unname(roi[[field]]),
          unname(expected_roi[[field]]),
          info = paste(mode, backend, field)
        )
      }
      expect_true(
        isTRUE(roi$unique_exclusive_fate_sum_matches_focal),
        info = paste(mode, backend, "fate partition")
      )
      expect_equal(roi$target_intersection_mode, mode)
      backend_results[[mode]][[backend]] <- list(result = result, roi = roi)
    }

    expect_equal(
      ftq_target_assignment_core(backend_results[[mode]]$r$result),
      ftq_target_assignment_core(backend_results[[mode]]$cpp$result),
      info = paste(mode, "complete per-voxel R/C++ parity")
    )
    r_roi <- backend_results[[mode]]$r$roi
    cpp_roi <- backend_results[[mode]]$cpp$roi
    r_roi$backend <- NULL
    cpp_roi$backend <- NULL
    expect_equal(
      r_roi,
      cpp_roi,
      info = paste(mode, "complete ROI-summary R/C++ parity")
    )
  }
})

test_that("endpoint voxels use the first focal-positive discretization", {
  dims <- c(8L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 6L)
  track <- list(rbind(c(2.6, 1, 1), c(6.2, 1, 1)))

  for (backend in c("r", "cpp")) {
    round_only <- ftq_target_assignment_run(
      focal,
      list(T = target),
      track,
      mode = "endpoint_voxels",
      backend = backend,
      discretization = "round"
    )
    any_method <- ftq_target_assignment_run(
      focal,
      list(T = target),
      track,
      mode = "endpoint_voxels",
      backend = backend,
      discretization = "any"
    )

    expect_equal(round_only$tracks, 0L, info = backend)
    expect_equal(round_only$T, 0L, info = backend)
    expect_equal(any_method$tracks, 1L, info = backend)
    expect_equal(any_method$T, 1L, info = backend)
    expect_equal(
      attr(any_method, "fibrtraq_discretization"),
      c("round", "floor", "ceil"),
      info = backend
    )
  }
})

test_that("local search is focal-conditioned and includes the physical radius boundary", {
  dims <- c(10L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, 3L)
  target <- ftq_target_assignment_mask(dims, 7L)
  tracks <- list(
    rbind(c(1, 1, 1), c(3, 1, 1), c(9, 1, 1)),
    rbind(c(1, 1, 1), c(9, 1, 1))
  )

  for (backend in c("r", "cpp")) {
    inclusive <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "local_search", backend = backend, radius_mm = 2
    )
    inside_only <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "local_search", backend = backend,
      radius_mm = 2 - 1e-7
    )

    # Only the streamline that intersects the focal voxel enters the counting
    # engine, even though both endpoint pairs are within 2 mm of the nodes.
    expect_equal(inclusive$tracks, 1L, info = backend)
    expect_equal(inclusive$T, 1L, info = backend)
    expect_equal(inside_only$tracks, 1L, info = backend)
    expect_equal(inside_only$T, 0L, info = backend)
  }
})

test_that("local search measures anisotropic scanner-space distance", {
  dims <- c(7L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 5L)
  tracks <- list(
    rbind(c(1.4, 1, 1), c(2, 1, 1), c(5.6, 1, 1)),
    rbind(c(1.3, 1, 1), c(2, 1, 1), c(5.7, 1, 1))
  )
  anisotropic_basis <- diag(c(3, 1, 1), nrow = 3L)

  for (backend in c("r", "cpp")) {
    anisotropic <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "local_search", backend = backend,
      voxel_to_world = anisotropic_basis
    )
    isotropic <- ftq_target_assignment_run(
      focal, list(T = target), tracks,
      mode = "local_search", backend = backend,
      voxel_to_world = diag(3)
    )

    expect_equal(anisotropic$tracks, 2L, info = backend)
    expect_equal(anisotropic$T, 1L, info = backend)
    expect_equal(isotropic$T, 2L, info = backend)
  }
})

test_that("NIfTI local search handles oblique geometry, flips, and off-grid endpoints", {
  skip_if_not_installed("RNifti")

  dims <- c(8L, 3L, 3L)
  focal <- array(0, dim = dims)
  focal[7, 2, 2] <- 1
  target <- array(0, dim = dims)
  target[1, 2, 2] <- 1
  track <- list(rbind(
    c(2, 2, 2),
    c(4, 2, 2),
    c(8.6, 2, 2)
  ))
  reference <- RNifti::asNifti(
    array(0, dim = dims),
    reference = list(
      sform_code = 2L,
      qform_code = 0L,
      srow_x = c(2, 0.4, 0.2, 0),
      srow_y = c(0, 1.5, 0.3, 0),
      srow_z = c(0, 0, 2.5, 0)
    )
  )
  RNifti::pixunits(reference) <- "mm"
  reference_path <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(reference, reference_path)

  geometry <- FibrTraQ:::ftq_local_search_geometry(
    focal,
    ref_image = reference_path
  )
  expect_true(any(abs(geometry$voxel_to_world[row(
    geometry$voxel_to_world
  ) != col(geometry$voxel_to_world)]) > 0))
  expect_equal(geometry$spatial_units, "mm")
  expect_equal(geometry$spatial_units_scale_to_mm, 1)

  backend_tables <- list()
  for (backend in c("r", "cpp")) {
    result <- FibrTraQ(
      fNode = focal,
      localROI = array(0, dim = dims),
      cNode_list = list(T = target),
      cNode_labels = "T",
      all_tracks = track,
      discretization = "round",
      track_space = "voxel",
      ref_image = reference_path,
      axis_flip = c(TRUE, FALSE, FALSE),
      assignment_mode = "local_search",
      search_radius_mm = 2,
      endpoint_analysis = "exact_and_dilated1",
      backend = backend
    )
    endpoint <- ftq_endpoint_tables(result)

    expect_equal(result$tracks, 1L, info = backend)
    expect_equal(result$T, 1L, info = backend)
    expect_equal(
      endpoint$summary$target_intersection_streamlines,
      c(0, 0),
      info = backend
    )
    expect_equal(
      attr(result, "fibrtraq_local_search_spatial_units"),
      "mm",
      info = backend
    )
    backend_tables[[backend]] <- endpoint
  }
  expect_equal(backend_tables$r, backend_tables$cpp, ignore_attr = TRUE)
})

test_that("local-search NIfTI spatial units are normalized to millimetres", {
  skip_if_not_installed("RNifti")

  dims <- c(12L, 2L, 2L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 10L)
  track <- list(rbind(c(2, 1, 1), c(7, 1, 1)))
  fixtures <- list(
    m = list(spacing = 0.001, scale = 1000, known = TRUE),
    mm = list(spacing = 1, scale = 1, known = TRUE),
    um = list(spacing = 1000, scale = 0.001, known = TRUE),
    Unknown = list(spacing = 1, scale = 1, known = FALSE)
  )

  for (unit in names(fixtures)) {
    fixture <- fixtures[[unit]]
    reference <- RNifti::asNifti(array(0, dim = dims))
    RNifti::pixdim(reference) <- rep(fixture$spacing, 3L)
    RNifti::pixunits(reference) <- unit
    reference_path <- tempfile(fileext = ".nii.gz")
    RNifti::writeNifti(reference, reference_path)

    geometry <- FibrTraQ:::ftq_local_search_geometry(
      focal,
      ref_image = reference_path
    )
    expect_equal(
      geometry$voxel_size_mm,
      c(1, 1, 1),
      tolerance = 1e-6,
      info = unit
    )
    expect_equal(
      geometry$spatial_units_scale_to_mm,
      fixture$scale,
      info = unit
    )
    expect_identical(
      geometry$spatial_units_known,
      fixture$known,
      info = unit
    )

    for (backend in c("r", "cpp")) {
      run <- function() {
        FibrTraQ(
          fNode = focal,
          localROI = array(0, dim = dims),
          cNode_list = list(T = target),
          cNode_labels = "T",
          all_tracks = track,
          discretization = "round",
          track_space = "voxel",
          ref_image = reference_path,
          axis_flip = c(FALSE, FALSE, FALSE),
          assignment_mode = "local_search",
          search_radius_mm = 2,
          endpoint_analysis = "none",
          backend = backend
        )
      }
      if (identical(unit, "Unknown")) {
        expect_warning(
          result <- run(),
          "does not declare NIfTI spatial units",
          fixed = TRUE,
          info = backend
        )
      } else {
        result <- run()
      }
      expect_equal(result$T, 0L, info = paste(unit, backend))
      expect_equal(
        attr(result, "fibrtraq_local_search_spatial_units_scale_to_mm"),
        fixture$scale,
        info = paste(unit, backend)
      )
    }
  }

  make_reference <- function(spacing, unit) {
    image <- RNifti::asNifti(array(0, dim = c(3L, 3L, 3L)))
    RNifti::pixdim(image) <- rep(spacing, 3L)
    RNifti::pixunits(image) <- unit
    image
  }
  mm_reference <- make_reference(2, "mm")
  equivalent_m_reference <- make_reference(0.002, "m")
  mismatched_m_reference <- make_reference(2, "m")
  expect_true(FibrTraQ:::ftq_local_search_affine_equal(
    mm_reference,
    equivalent_m_reference
  ))
  expect_false(FibrTraQ:::ftq_local_search_affine_equal(
    mm_reference,
    mismatched_m_reference
  ))
  expect_error(
    FibrTraQ:::ftq_local_search_geometry(
      mm_reference,
      ref_image = mismatched_m_reference
    ),
    "must use the same affine/grid",
    fixed = TRUE
  )
})

test_that("world and TrackVis millimetres map through non-mm NIfTI units", {
  skip_if_not_installed("RNifti")

  dims <- c(12L, 2L, 2L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 10L)
  reference <- RNifti::asNifti(array(0, dim = dims))
  RNifti::pixdim(reference) <- rep(0.001, 3L)
  RNifti::pixunits(reference) <- "m"
  reference_path <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(reference, reference_path)
  endpoints_mm <- FibrTraQ:::ftq_voxel_to_world_mm(
    rbind(c(2, 1, 1), c(10, 1, 1)),
    reference_path
  )

  for (backend in c("r", "cpp")) {
    for (space in c("world", "trackvis")) {
      args <- list(
        fNode = focal,
        localROI = array(0, dim = dims),
        cNode_list = list(T = target),
        cNode_labels = "T",
        all_tracks = list(endpoints_mm),
        discretization = "round",
        track_space = space,
        ref_image = reference_path,
        axis_flip = c(FALSE, FALSE, FALSE),
        assignment_mode = "local_search",
        search_radius_mm = 2,
        endpoint_analysis = "none",
        backend = backend
      )
      if (identical(space, "trackvis")) args$track_affine <- diag(4)
      result <- do.call(FibrTraQ, args)

      expect_equal(result$tracks, 1L, info = paste(space, backend))
      expect_equal(result$T, 1L, info = paste(space, backend))
      expect_equal(
        attr(result, "fibrtraq_track_space"),
        space,
        info = paste(space, backend)
      )
    }
  }
})

test_that("local search resolves nearest labels and rejects ambiguous ties and overlaps", {
  dims <- c(10L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target_a <- ftq_target_assignment_mask(dims, 7L)
  target_b <- ftq_target_assignment_mask(dims, 8L)
  nearest_track <- list(rbind(c(2, 1, 1), c(4, 1, 1), c(7.4, 1, 1)))
  tie_track <- list(rbind(c(2, 1, 1), c(4, 1, 1), c(7.5, 1, 1)))

  for (backend in c("r", "cpp")) {
    nearest <- ftq_target_assignment_run(
      focal, list(A = target_a, B = target_b), nearest_track,
      mode = "local_search", backend = backend
    )
    nearest_reordered <- ftq_target_assignment_run(
      focal, list(B = target_b, A = target_a), nearest_track,
      mode = "local_search", backend = backend
    )
    tied <- ftq_target_assignment_run(
      focal, list(A = target_a, B = target_b), tie_track,
      mode = "local_search", backend = backend
    )
    same_label_tie <- ftq_target_assignment_run(
      focal,
      list(A = ftq_target_assignment_mask(dims, 7:8)),
      tie_track,
      mode = "local_search",
      backend = backend
    )
    resolved_near_tie <- ftq_target_assignment_run(
      focal,
      list(A = target_a, B = target_b),
      list(rbind(c(2, 1, 1), c(4, 1, 1), c(7.500000000001, 1, 1))),
      mode = "local_search",
      backend = backend
    )
    wide_radius <- ftq_target_assignment_run(
      focal,
      list(A = target_a, B = target_b),
      list(rbind(c(2, 1, 1), c(4, 1, 1), c(7.49, 1, 1))),
      mode = "local_search",
      backend = backend,
      radius_mm = 1e6
    )

    expect_equal(c(nearest$A, nearest$B), c(1L, 0L), info = backend)
    expect_equal(
      c(A = nearest_reordered$A, B = nearest_reordered$B),
      c(A = 1L, B = 0L),
      info = backend
    )
    expect_equal(c(tied$A, tied$B), c(0L, 0L), info = backend)
    expect_equal(same_label_tie$A, 1L, info = backend)
    expect_equal(
      c(resolved_near_tie$A, resolved_near_tie$B),
      c(0L, 1L),
      info = backend
    )
    expect_equal(c(wide_radius$A, wide_radius$B), c(1L, 0L),
                 info = backend)

    overlap_a <- ftq_target_assignment_mask(dims, 7L)
    overlap_b <- ftq_target_assignment_mask(dims, 7L)
    overlap_track <- list(rbind(c(2, 1, 1), c(7, 1, 1)))
    overlap_full <- ftq_target_assignment_run(
      focal, list(A = overlap_a, B = overlap_b), overlap_track,
      mode = "full_streamline", backend = backend
    )
    overlap_endpoint <- ftq_target_assignment_run(
      focal, list(A = overlap_a, B = overlap_b), overlap_track,
      mode = "endpoint_voxels", backend = backend
    )
    overlap_local <- ftq_target_assignment_run(
      focal, list(A = overlap_a, B = overlap_b), overlap_track,
      mode = "local_search", backend = backend
    )
    expect_equal(c(overlap_full$A, overlap_full$B), c(1L, 1L),
                 info = backend)
    expect_equal(c(overlap_endpoint$A, overlap_endpoint$B), c(0L, 0L),
                 info = backend)
    expect_equal(c(overlap_local$A, overlap_local$B), c(0L, 0L),
                 info = backend)
  }
})

test_that("local search does not stop at the first labelled coordinate shell", {
  anisotropic_dims <- c(30L, 4L, 2L)
  anisotropic_basis <- diag(c(0.1, 5, 2))
  anisotropic_focal <- ftq_target_assignment_point_mask(
    anisotropic_dims,
    rbind(c(30L, 4L, 1L))
  )
  anisotropic_inner <- ftq_target_assignment_point_mask(
    anisotropic_dims,
    rbind(c(2L, 3L, 1L))
  )
  anisotropic_outer <- ftq_target_assignment_point_mask(
    anisotropic_dims,
    rbind(c(22L, 2L, 1L))
  )
  anisotropic_track <- list(rbind(c(30, 4, 1), c(2, 2, 1)))

  shear_dims <- c(30L, 30L, 3L)
  shear_basis <- rbind(
    c(1, 0, 0),
    c(0.99, 0.01, 0),
    c(0, 0, 1)
  )
  shear_focal <- ftq_target_assignment_point_mask(
    shear_dims,
    rbind(c(30L, 30L, 3L))
  )
  shear_inner <- ftq_target_assignment_point_mask(
    shear_dims,
    rbind(c(16L, 15L, 2L))
  )
  shear_outer <- ftq_target_assignment_point_mask(
    shear_dims,
    rbind(c(25L, 5L, 2L))
  )
  shear_track <- list(rbind(c(30, 30, 3), c(15, 15, 2)))

  tie_dims <- c(25L, 5L, 2L)
  tie_basis <- diag(c(0.1, 1, 1))
  tie_focal <- ftq_target_assignment_point_mask(
    tie_dims,
    rbind(c(25L, 5L, 2L))
  )
  tie_inner <- ftq_target_assignment_point_mask(
    tie_dims,
    rbind(c(5L, 3L, 1L))
  )
  tie_outer <- ftq_target_assignment_point_mask(
    tie_dims,
    rbind(c(15L, 2L, 1L))
  )
  tie_track <- list(rbind(c(25, 5, 2), c(5, 2, 1)))

  for (backend in c("r", "cpp")) {
    anisotropic <- ftq_target_assignment_run(
      anisotropic_focal,
      list(inner = anisotropic_inner, outer = anisotropic_outer),
      anisotropic_track,
      mode = "local_search",
      backend = backend,
      radius_mm = 6,
      voxel_to_world = anisotropic_basis
    )
    expect_equal(
      c(inner = anisotropic$inner, outer = anisotropic$outer),
      c(inner = 0L, outer = 1L),
      info = paste(backend, "anisotropic outer-nearest")
    )

    sheared <- ftq_target_assignment_run(
      shear_focal,
      list(inner = shear_inner, outer = shear_outer),
      shear_track,
      mode = "local_search",
      backend = backend,
      radius_mm = 2,
      voxel_to_world = shear_basis
    )
    expect_equal(
      c(inner = sheared$inner, outer = sheared$outer),
      c(inner = 0L, outer = 1L),
      info = paste(backend, "shear cancellation")
    )

    tied <- ftq_target_assignment_run(
      tie_focal,
      list(inner = tie_inner, outer = tie_outer),
      tie_track,
      mode = "local_search",
      backend = backend,
      radius_mm = 2,
      voxel_to_world = tie_basis
    )
    expect_equal(
      c(inner = tied$inner, outer = tied$outer),
      c(inner = 0L, outer = 0L),
      info = paste(backend, "cross-shell physical tie")
    )
  }
})

test_that("endpoint and local assignments are reversal and voxelization invariant", {
  dims <- c(9L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 7L)
  mode_tracks <- list(
    endpoint_voxels = list(rbind(c(2, 1, 1), c(4, 1, 1), c(7, 1, 1))),
    local_search = list(
      rbind(c(1, 1, 1), c(2, 1, 1), c(4, 1, 1), c(8, 1, 1))
    )
  )

  for (backend in c("r", "cpp")) {
    for (mode in names(mode_tracks)) {
      baseline <- ftq_target_assignment_run(
        focal, list(T = target), mode_tracks[[mode]],
        mode = mode, backend = backend
      )
      reversed <- ftq_target_assignment_run(
        focal, list(T = target), mode_tracks[[mode]],
        mode = mode, backend = backend, reverse_streamlines = TRUE
      )
      interpolated <- ftq_target_assignment_run(
        focal, list(T = target), mode_tracks[[mode]],
        mode = mode, backend = backend,
        voxelization = "segment_interpolation"
      )

      expect_equal(baseline$T, 1L, info = paste(backend, mode))
      expect_equal(
        ftq_target_assignment_core(reversed),
        ftq_target_assignment_core(baseline),
        info = paste(backend, mode, "reversal")
      )
      expect_equal(
        ftq_target_assignment_core(interpolated),
        ftq_target_assignment_core(baseline),
        info = paste(backend, mode, "voxelization")
      )
    }
  }
})

test_that("endpoint audit keeps a full-path denominator under endpoint assignment", {
  dims <- c(9L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 6L)
  track <- list(rbind(c(2, 1, 1), c(6, 1, 1), c(9, 1, 1)))
  backend_tables <- list()

  for (backend in c("r", "cpp")) {
    result <- ftq_target_assignment_run(
      focal, list(T = target), track,
      mode = "endpoint_voxels", backend = backend,
      endpoint_analysis = "exact_and_dilated1"
    )
    roi <- attr(result, "fibrtraq_roi_summary")
    endpoint <- ftq_endpoint_tables(result)
    checks <- FibrTraQ:::ftq_validate_endpoint_tables(
      endpoint$summary,
      endpoint$penetration,
      endpoint$context,
      roi_summary = roi
    )$checks
    denominator_check <- checks[
      checks$check == "endpoint_denominator_preserved_as_full_path",
      ,
      drop = FALSE
    ]

    expect_equal(result$T, 0L, info = backend)
    expect_equal(roi$`unique_target_streamlines:T`, 0, info = backend)
    expect_equal(
      roi$`unique_full_path_target_streamlines:T`,
      1,
      info = backend
    )
    expect_equal(
      endpoint$summary$target_intersection_streamlines,
      c(1, 1),
      info = backend
    )
    expect_true(denominator_check$passed, info = backend)
    expect_match(
      denominator_check$detail,
      "reconciled to independently persisted full-streamline counts",
      fixed = TRUE
    )

    tampered_roi <- roi
    tampered_roi$`unique_full_path_target_streamlines:T` <- 0
    tampered_checks <- FibrTraQ:::ftq_validate_endpoint_tables(
      endpoint$summary,
      endpoint$penetration,
      endpoint$context,
      roi_summary = tampered_roi
    )$checks
    tampered_denominator_check <- tampered_checks[
      tampered_checks$check ==
        "endpoint_denominator_preserved_as_full_path",
      ,
      drop = FALSE
    ]
    expect_false(tampered_denominator_check$passed, info = backend)
    backend_tables[[backend]] <- endpoint
  }
  expect_equal(
    backend_tables$r,
    backend_tables$cpp,
    ignore_attr = TRUE
  )
})

test_that("invalid target modes, radii, and missing local geometry are rejected", {
  dims <- c(5L, 1L, 1L)
  focal <- ftq_target_assignment_mask(dims, 2L)
  target <- ftq_target_assignment_mask(dims, 4L)
  track <- list(rbind(c(2, 1, 1), c(4, 1, 1)))

  expect_error(
    ftq_target_assignment_run(
      focal, list(T = target), track, mode = "not_a_mode"
    ),
    "target_intersection_mode must be one of",
    fixed = TRUE
  )
  for (radius in list(0, -1, Inf, NA_real_, c(1, 2))) {
    expect_error(
      ftq_target_assignment_run(
        focal, list(T = target), track,
        mode = "local_search", radius_mm = radius
      ),
      "local_search_radius_mm must be one finite number greater than zero",
      fixed = TRUE
    )
  }

  expect_error(
    FibrTraQ(
      fNode = focal,
      localROI = array(0, dim = dims),
      cNode_list = list(T = target),
      cNode_labels = "T",
      all_tracks = track,
      discretization = "round",
      track_space = "voxel",
      axis_flip = c(FALSE, FALSE, FALSE),
      assignment_mode = "local_search",
      endpoint_analysis = "none",
      backend = "r"
    ),
    "requires reference-image spatial metadata",
    fixed = TRUE
  )

  bad_reference <- RNifti::asNifti(array(0, dim = c(6L, 2L, 2L)))
  RNifti::pixunits(bad_reference) <- "mm"
  bad_reference_path <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(bad_reference, bad_reference_path)
  expect_error(
    FibrTraQ(
      fNode = focal,
      localROI = array(0, dim = dims),
      cNode_list = list(T = target),
      cNode_labels = "T",
      all_tracks = track,
      discretization = "round",
      track_space = "voxel",
      ref_image = bad_reference_path,
      axis_flip = c(FALSE, FALSE, FALSE),
      assignment_mode = "local_search",
      endpoint_analysis = "none",
      backend = "r"
    ),
    "ref_image dimensions must match the fNode grid",
    fixed = TRUE
  )

  expect_error(
    ftq_target_assignment_run(
      focal,
      list(T = target),
      track,
      mode = "local_search",
      radius_mm = 1e308,
      voxel_to_world = diag(c(0.1, 1, 1))
    ),
    "produces a non-finite voxel search extent",
    fixed = TRUE
  )
})

test_that("manifests and file wrapper propagate target assignment settings", {
  root <- file.path(tempdir(), "target-assignment-manifest")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  manifest <- list(
    dataset_root = root,
    subject = "sub-test",
    description = "target-assignment",
    fnode = "missing-focal.nii.gz",
    tractogram = "missing-tractogram.trk",
    targets = list(T = "missing-target.nii.gz"),
    analysis = list(
      target_intersection_mode = "end_voxels",
      local_search_radius_mm = 3.25,
      endpoint_analysis = "none"
    ),
    outputs = list(
      output_dir = file.path(root, "output"),
      write_maps = FALSE,
      write_report = FALSE,
      overwrite = TRUE
    )
  )
  parsed <- ftq_read_manifest(manifest)
  expect_equal(parsed$target_intersection_mode, "endpoint_voxels")
  expect_equal(parsed$local_search_radius_mm, 3.25)

  preflight <- ftq_preflight(parsed)
  expect_true(preflight$checks$passed[
    preflight$checks$check == "target_intersection_mode_supported"
  ])
  expect_true(preflight$checks$passed[
    preflight$checks$check == "local_search_radius_mm_positive"
  ])

  invalid <- manifest
  invalid$analysis$local_search_radius_mm <- 0
  expect_error(
    ftq_read_manifest(invalid),
    "local_search_radius_mm must be one finite number greater than zero",
    fixed = TRUE
  )

  inputs <- file.path(root, c("focal.nii.gz", "target.nii.gz", "tracks.trk"))
  invisible(file.create(inputs))
  expect_identical(
    names(formals(ftq_run_files))[pmatch("target", names(formals(ftq_run_files)))],
    "target_masks"
  )
  expect_identical(
    names(formals(ftq_run_files))[pmatch("local", names(formals(ftq_run_files)))],
    "local_mask"
  )
  expect_identical(
    names(formals(ftq_run_regional_study))[
      pmatch("local", names(formals(ftq_run_regional_study)))
    ],
    "local_masks"
  )
  captured <- NULL
  local_mocked_bindings(
    ftq_run_regional_study = function(...) {
      captured <<- list(...)
      structure(list(captured = captured), class = "ftq_regional_study_result")
    },
    .package = "FibrTraQ"
  )
  wrapper_result <- ftq_run_files(
    tractogram = inputs[[3L]],
    focal_mask = inputs[[1L]],
    target_masks = c(T = inputs[[2L]]),
    output_dir = file.path(root, "wrapper-output"),
    write_maps = FALSE,
    write_report = FALSE,
    write_exploratory_report = FALSE,
    write_qc_viewer = FALSE,
    verbose = FALSE,
    progress = FALSE,
    assignment_mode = "local_search",
    search_radius_mm = 3.25
  )
  expect_s3_class(wrapper_result, "ftq_regional_study_result")
  expect_equal(captured$assignment_mode, "local_search")
  expect_equal(captured$search_radius_mm, 3.25)
})
