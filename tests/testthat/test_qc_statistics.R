library(testthat)
library(FibrTraQ)

ftq_write_qc_test_trk <- function(path, tracks, dims = c(5L, 5L, 5L)) {
  write_fixed <- function(con, text, n) {
    bytes <- rep(as.raw(0), n)
    text_bytes <- charToRaw(text)
    keep <- seq_len(min(length(text_bytes), n))
    bytes[keep] <- text_bytes[keep]
    writeBin(bytes, con)
  }
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
  writeBin(as.numeric(c(1, 0, 0, 0, 1, 0)), con, size = 4L, endian = endian)
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

ftq_scene_payload <- function(path) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  m <- regexpr("<script id=\"scene-data\" type=\"application/json\">.*?</script>", html)
  stopifnot(m[1] > 0)
  block <- regmatches(html, m)
  json <- sub("^<script id=\"scene-data\" type=\"application/json\">", "", block)
  json <- sub("</script>$", "", json)
  list(html = html, payload = jsonlite::fromJSON(json, simplifyVector = FALSE))
}

test_that("interactive QC viewer renders full focal volume and segment controls", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "interactive-qc-viewer")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  fnode <- file.path(tmp, "fnode.nii.gz")
  target_a <- file.path(tmp, "target-a.nii.gz")
  target_b <- file.path(tmp, "target-b.nii.gz")
  background <- file.path(tmp, "background.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")

  fnode_arr <- array(0, dim = c(5, 5, 5))
  fnode_arr[2, 2, 2] <- 1
  fnode_arr[2, 2, 3] <- 1
  a_arr <- array(0, dim = c(5, 5, 5))
  a_arr[5, 2, 2] <- 1
  b_arr <- array(0, dim = c(5, 5, 5))
  b_arr[5, 2, 3] <- 1
  RNifti::writeNifti(fnode_arr, fnode)
  RNifti::writeNifti(a_arr, target_a)
  RNifti::writeNifti(b_arr, target_b)
  RNifti::writeNifti(array(seq_len(125), dim = c(5, 5, 5)), background)
  ftq_write_qc_test_trk(
    tractogram,
    list(
      rbind(c(1, 2, 2), c(2, 2, 2), c(5, 2, 2)),
      rbind(c(1, 2, 3), c(2, 2, 3), c(5, 2, 3)),
      rbind(c(1, 1, 1), c(1, 1, 2))
    )
  )

  manifest <- list(
    dataset_root = tmp,
    subject = "sub-test",
    description = "qc",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(A = basename(target_a), B = basename(target_b)),
    background_image = basename(background),
    coordinate = list(track_space = "voxel", discretization = "round", axis_flip = c(FALSE, FALSE, FALSE)),
    execution = list(track_reader = "stream", track_chunk_size = 2L),
    outputs = list(output_dir = file.path(tmp, "out"), write_maps = FALSE, write_report = FALSE, overwrite = TRUE)
  )
  viewer <- ftq_interactive_qc_viewer(
    manifest,
    output_html = file.path(tmp, "viewer.html"),
    sample_per_target = 1L,
    overwrite = TRUE
  )
  expect_true(file.exists(viewer$html))
  expect_equal(viewer$focal_voxels, 2L)
  scene <- ftq_scene_payload(viewer$html)
  expect_equal(scene$payload$focalCount, 2L)
  expect_equal(length(scene$payload$focal), 2L)
  expect_equal(length(scene$payload$targets[[1]]$boundary), 1L)
  expect_equal(scene$payload$targets[[1]]$voxelCount, 1L)
  expect_equal(length(scene$payload$streamlines), 2L)
  expect_true(length(scene$payload$streamlines[[1]]$hits) >= 1L)
  expect_true(length(scene$payload$streamlines[[1]]$targetHits) >= 1L)
  expect_true(length(scene$payload$streamlines[[1]]$targetSegment) >= 1L)
  expect_true(scene$payload$endpointAnalysis)
  expect_true(scene$payload$streamlines[[1]]$endpointClassExact %in% scene$payload$endpointClasses)
  expect_true(scene$payload$streamlines[[1]]$endpointClassDilated1 %in% scene$payload$endpointClasses)
  expect_true(length(scene$payload$streamlines[[1]]$endpoints) >= 1L)
  expect_true(grepl('data-mode="focal"', scene$html, fixed = TRUE))
  expect_true(grepl('data-mode="targetSegment"', scene$html, fixed = TRUE))
  expect_true(grepl('data-mode="after"', scene$html, fixed = TRUE))
  expect_true(grepl("cNode Extents", scene$html, fixed = TRUE))
  expect_true(grepl("cNode-Hit Markers", scene$html, fixed = TRUE))
  expect_true(grepl("function hits", scene$html, fixed = TRUE))
  expect_true(grepl("function targetHits", scene$html, fixed = TRUE))
  expect_true(grepl("1-voxel sensitivity", scene$html, fixed = TRUE))
  expect_true(grepl("classFilters", scene$html, fixed = TRUE))
  expect_true(grepl("function endpointMarkers", scene$html, fixed = TRUE))
  expect_identical(scene$payload$samplingStrategy, "target_quota_scan_order")
  expect_match(scene$payload$metadata$subtitle, "Target-balanced quota QC subset",
               fixed = TRUE)
  expect_match(scene$payload$metadata$subtitle, "not a random or prevalence-weighted sample",
               fixed = TRUE)
  expect_false(grepl("representative QC sample", scene$html, fixed = TRUE))
  expect_true(grepl("Advanced controls", scene$html, fixed = TRUE))
  expect_true(grepl(
    '<button id="toggleHits" type="button" aria-pressed="false">',
    scene$html,
    fixed = TRUE
  ))
  expect_true(grepl(
    '<button id="toggleTargetHits" type="button" aria-pressed="false">',
    scene$html,
    fixed = TRUE
  ))
  expect_true(grepl(
    "let x=pt[0]-c[0],y=pt[2]-c[2],z=-(pt[1]-c[1])",
    scene$html,
    fixed = TRUE
  ))
  expect_true(grepl('[["R",', scene$html, fixed = TRUE))
  expect_true(grepl('["S",', scene$html, fixed = TRUE))

  quota_viewer <- ftq_interactive_qc_viewer(
    manifest,
    output_html = file.path(tmp, "viewer-quota.html"),
    target_labels = c("A", "B"),
    sample_per_target = c(A = 1L, B = 0L),
    overwrite = TRUE
  )
  quota_scene <- ftq_scene_payload(quota_viewer$html)
  expect_equal(length(quota_scene$payload$streamlines), 1L)
  expect_equal(quota_scene$payload$streamlines[[1]]$target, "A")
  expect_equal(
    quota_viewer$summary$sampled_streamlines[quota_viewer$summary$target == "B"],
    0L
  )

  legacy_manifest <- ftq_read_manifest(manifest)
  legacy_manifest$target_intersection_mode <- NULL
  legacy_manifest$local_search_radius_mm <- NULL
  legacy_viewer <- ftq_interactive_qc_viewer(
    list(
      manifest = legacy_manifest,
      outputs = list(context = NULL),
      result = NULL
    ),
    output_html = file.path(tmp, "viewer-legacy-run.html"),
    target_labels = "A",
    sample_per_target = 1L,
    overwrite = TRUE
  )
  expect_equal(legacy_viewer$target_intersection_mode, "full_streamline")
  expect_true(file.exists(legacy_viewer$html))
})

test_that("interactive QC viewer treats target labels as data, not markup", {
  malicious_label <- paste0(
    "</script><script>globalThis.FIBRTRAQXSS=1</script>",
    "<img src=x onerror=alert(1)>"
  )
  output_html <- tempfile("fibrtraq-qc-untrusted-label-", fileext = ".html")
  points <- list(c(0, 0, 0), c(1, 1, 1))
  streamlines <- list(list(
    target = malicious_label,
    points = points,
    focal = points,
    target_segment = points,
    after = points,
    hits = list(c(0, 0, 0)),
    target_hits = list(c(1, 1, 1)),
    endpoint_annotation = NULL
  ))
  target_centers <- data.frame(
    target = malicious_label,
    x = 1,
    y = 1,
    z = 1,
    stringsAsFactors = FALSE
  )
  target_boundaries <- stats::setNames(list(list(
    points = data.frame(x = 1, y = 1, z = 1),
    voxel_count = 1L,
    rendered_count = 1L,
    render_mode = "full"
  )), malicious_label)
  colors <- stats::setNames("#112233", malicious_label)
  sample_counts <- stats::setNames(1L, malicious_label)

  FibrTraQ:::ftq_qc_write_html(
    path = output_html,
    streamlines = streamlines,
    focal_world = matrix(c(0, 0, 0), nrow = 1L),
    focal_total_count = 1L,
    focal_rendered_count = 1L,
    focal_render_mode = "full",
    target_centers = target_centers,
    target_boundaries = target_boundaries,
    background_world = data.frame(x = 0, y = 0, z = 0),
    context_label = "none",
    colors = colors,
    scanned = 1L,
    focal_hits = 1L,
    sample_counts = sample_counts,
    endpoint_analysis = FALSE,
    title = "Untrusted target-label regression"
  )

  scene <- ftq_scene_payload(output_html)
  expect_identical(scene$payload$legend[[1L]]$target, malicious_label)
  expect_false(grepl(malicious_label, scene$html, fixed = TRUE))
  escaped_script_end <- paste0(
    intToUtf8(92L),
    "u003c",
    intToUtf8(92L),
    "/script",
    intToUtf8(92L),
    "u003e"
  )
  expect_true(grepl(escaped_script_end, scene$html, fixed = TRUE))
  expect_false(grepl("elLegend.innerHTML", scene$html, fixed = TRUE))
  expect_true(grepl("document.createTextNode(String(d.target))", scene$html,
                    fixed = TRUE))
})

test_that("QC target subsets retain omitted labels for assignment ambiguity", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "interactive-qc-target-subset-assignment")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  fnode <- file.path(tmp, "fnode.nii.gz")
  target_a <- file.path(tmp, "target-a.nii.gz")
  target_b <- file.path(tmp, "target-b.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")

  fnode_arr <- array(0, dim = c(6, 3, 3))
  fnode_arr[2, 2, 2] <- 1
  target_arr <- array(0, dim = c(6, 3, 3))
  target_arr[6, 2, 2] <- 1
  RNifti::writeNifti(fnode_arr, fnode)
  RNifti::writeNifti(target_arr, target_a)
  RNifti::writeNifti(target_arr, target_b)
  ftq_write_qc_test_trk(
    tractogram,
    list(rbind(c(2, 2, 2), c(5, 2, 2))),
    dims = c(6L, 3L, 3L)
  )

  base_manifest <- list(
    dataset_root = tmp,
    subject = "sub-test",
    description = "qc-target-subset",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(A = basename(target_a), B = basename(target_b)),
    coordinate = list(
      track_space = "voxel",
      discretization = "round",
      axis_flip = c(FALSE, FALSE, FALSE)
    ),
    analysis = list(
      target_intersection_mode = "local_search",
      local_search_radius_mm = 2,
      endpoint_analysis = "none"
    ),
    outputs = list(
      output_dir = file.path(tmp, "out"),
      write_maps = FALSE,
      write_report = FALSE,
      overwrite = TRUE
    )
  )

  for (reader in c("memory", "stream")) {
    manifest <- base_manifest
    manifest$execution <- list(
      track_reader = reader,
      track_chunk_size = 1L
    )
    expect_error(
      ftq_interactive_qc_viewer(
        manifest,
        output_html = file.path(tmp, paste0("viewer-", reader, ".html")),
        target_labels = "A",
        sample_per_target = 1L,
        overwrite = TRUE
      ),
      "No sampled focal-positive streamlines were connected",
      fixed = TRUE,
      info = reader
    )
  }
})

test_that("QC sampling follows each selected target-assignment mode", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "interactive-qc-assignment-modes")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  fnode <- file.path(tmp, "fnode.nii.gz")
  target <- file.path(tmp, "target.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")
  endpoint_tractogram <- file.path(tmp, "endpoint-tracts.trk")

  fnode_arr <- array(0, dim = c(8, 3, 3))
  fnode_arr[2, 2, 2] <- 1
  target_arr <- array(0, dim = c(8, 3, 3))
  target_arr[7, 2, 2] <- 1
  fnode_nifti <- RNifti::asNifti(fnode_arr)
  target_nifti <- RNifti::asNifti(target_arr)
  RNifti::pixunits(fnode_nifti) <- "mm"
  RNifti::pixunits(target_nifti) <- "mm"
  RNifti::writeNifti(fnode_nifti, fnode)
  RNifti::writeNifti(target_nifti, target)
  ftq_write_qc_test_trk(
    tractogram,
    list(rbind(c(2, 2, 2), c(7, 2, 2), c(8, 2, 2))),
    dims = c(8L, 3L, 3L)
  )
  ftq_write_qc_test_trk(
    endpoint_tractogram,
    list(rbind(c(2, 2, 2), c(7, 2, 2))),
    dims = c(8L, 3L, 3L)
  )

  focal_image <- RNifti::readNifti(fnode, internal = FALSE)
  target_image <- RNifti::readNifti(target, internal = FALSE)
  interior_track <- list(rbind(c(2, 2, 2), c(7, 2, 2), c(8, 2, 2)))
  direct_counts <- vapply(
    c("full_streamline", "endpoint_voxels", "local_search"),
    function(mode) {
      result <- FibrTraQ(
        fNode = focal_image,
        localROI = focal_image,
        cNode_list = list(T = target_image),
        cNode_labels = "T",
        all_tracks = interior_track,
        discretization = "round",
        track_space = "voxel",
        axis_flip = c(FALSE, FALSE, FALSE),
        backend = "r",
        endpoint_analysis = "none",
        assignment_mode = mode,
        search_radius_mm = 1.1
      )
      unname(attr(result, "fibrtraq_roi_summary")[[
        "unique_target_streamlines:T"
      ]])
    },
    numeric(1)
  )
  expect_equal(unname(direct_counts), c(1, 0, 1))
  expected_target_world <- as.numeric(RNifti::voxelToWorld(
    matrix(c(7, 2, 2), ncol = 3),
    focal_image,
    simple = FALSE
  ))

  base_manifest <- list(
    dataset_root = tmp,
    subject = "sub-test",
    description = "qc-assignment-modes",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(T = basename(target)),
    coordinate = list(
      track_space = "voxel",
      discretization = "round",
      axis_flip = c(FALSE, FALSE, FALSE)
    ),
    outputs = list(
      output_dir = file.path(tmp, "out"),
      write_maps = FALSE,
      write_report = FALSE,
      overwrite = TRUE
    )
  )

  for (reader in c("memory", "stream")) {
    for (mode in c("full_streamline", "local_search")) {
      manifest <- base_manifest
      manifest$analysis <- list(
        target_intersection_mode = mode,
        local_search_radius_mm = 1.1,
        endpoint_analysis = "none"
      )
      manifest$execution <- list(
        track_reader = reader,
        track_chunk_size = 1L
      )
      viewer <- ftq_interactive_qc_viewer(
        manifest,
        output_html = file.path(
          tmp,
          paste0("viewer-", reader, "-", mode, ".html")
        ),
        sample_per_target = 1L,
        overwrite = TRUE
      )
      scene <- ftq_scene_payload(viewer$html)

      expect_equal(viewer$target_intersection_mode, mode,
                   info = paste(reader, mode))
      expect_equal(length(scene$payload$streamlines), 1L,
                   info = paste(reader, mode))
      expect_equal(scene$payload$streamlines[[1]]$target, "T",
                   info = paste(reader, mode))
      expect_equal(scene$payload$targetIntersectionMode, mode,
                   info = paste(reader, mode))
      expect_true(length(scene$payload$streamlines[[1]]$targetHits) >= 1L,
                  info = paste(reader, mode))
      expect_equal(
        unlist(
          scene$payload$streamlines[[1]]$targetHits[[1]],
          use.names = FALSE
        ),
        expected_target_world,
        tolerance = 1e-7,
        info = paste(reader, mode, "target marker")
      )
      if (identical(mode, "local_search")) {
        expect_equal(viewer$local_search_radius_mm, 1.1,
                     info = reader)
        expect_equal(scene$payload$localSearchRadiusMm, 1.1,
                     info = reader)
      } else {
        expect_equal(
          length(scene$payload$localSearchRadiusMm),
          0L,
          info = reader
        )
      }
    }

    endpoint_manifest <- base_manifest
    endpoint_manifest$analysis <- list(
      target_intersection_mode = "endpoint_voxels",
      local_search_radius_mm = 1.1,
      endpoint_analysis = "none"
    )
    endpoint_manifest$execution <- list(
      track_reader = reader,
      track_chunk_size = 1L
    )
    expect_error(
      ftq_interactive_qc_viewer(
        endpoint_manifest,
        output_html = file.path(
          tmp,
          paste0("viewer-", reader, "-endpoint-voxels.html")
        ),
        sample_per_target = 1L,
        overwrite = TRUE
      ),
      "No sampled focal-positive streamlines were connected",
      fixed = TRUE,
      info = reader
    )

    endpoint_positive_manifest <- endpoint_manifest
    endpoint_positive_manifest$tractogram <- basename(endpoint_tractogram)
    endpoint_viewer <- ftq_interactive_qc_viewer(
      endpoint_positive_manifest,
      output_html = file.path(
        tmp,
        paste0("viewer-", reader, "-endpoint-positive.html")
      ),
      sample_per_target = 1L,
      overwrite = TRUE
    )
    endpoint_scene <- ftq_scene_payload(endpoint_viewer$html)
    expect_equal(endpoint_viewer$target_intersection_mode, "endpoint_voxels",
                 info = reader)
    expect_equal(length(endpoint_scene$payload$streamlines), 1L,
                 info = reader)
    expect_equal(endpoint_scene$payload$streamlines[[1]]$target, "T",
                 info = reader)
    expect_equal(
      unlist(
        endpoint_scene$payload$streamlines[[1]]$targetHits[[1]],
        use.names = FALSE
      ),
      expected_target_world,
      tolerance = 1e-7,
      info = paste(reader, "endpoint target marker")
    )
  }
})

test_that("interactive QC viewer renders axis-flipped tracks in ROI image space", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "interactive-qc-viewer-axis-flip")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  fnode <- file.path(tmp, "fnode.nii.gz")
  target <- file.path(tmp, "target.nii.gz")
  local <- file.path(tmp, "local.nii.gz")
  background <- file.path(tmp, "background.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")

  fnode_arr <- array(0, dim = c(5, 5, 5))
  fnode_arr[4, 2, 2] <- 1
  target_arr <- array(0, dim = c(5, 5, 5))
  target_arr[1, 2, 2] <- 1
  local_arr <- array(0, dim = c(5, 5, 5))
  local_arr[4, 2, 2] <- 1
  local_arr[1, 2, 2] <- 1
  RNifti::writeNifti(fnode_arr, fnode)
  RNifti::writeNifti(target_arr, target)
  RNifti::writeNifti(local_arr, local)
  RNifti::writeNifti(array(1, dim = c(5, 5, 5)), background)
  ftq_write_qc_test_trk(
    tractogram,
    list(rbind(c(2, 2, 2), c(5, 2, 2)))
  )

  manifest <- list(
    dataset_root = tmp,
    subject = "sub-test",
    description = "qc-axis-flip",
    fnode = basename(fnode),
    local = basename(local),
    tractogram = basename(tractogram),
    targets = list(A = basename(target)),
    background_image = basename(background),
    coordinate = list(track_space = "voxel", discretization = "round", axis_flip = c(TRUE, FALSE, FALSE)),
    execution = list(track_reader = "stream", track_chunk_size = 2L),
    outputs = list(output_dir = file.path(tmp, "out"), write_maps = FALSE, write_report = FALSE, overwrite = TRUE)
  )
  viewer <- ftq_interactive_qc_viewer(
    manifest,
    output_html = file.path(tmp, "viewer.html"),
    sample_per_target = 1L,
    overwrite = TRUE
  )
  scene <- ftq_scene_payload(viewer$html)
  img <- RNifti::readNifti(fnode, internal = FALSE)
  expected_focal <- as.numeric(RNifti::voxelToWorld(matrix(c(4, 2, 2), ncol = 3), img, simple = FALSE))
  expected_target <- as.numeric(RNifti::voxelToWorld(matrix(c(1, 2, 2), ncol = 3), img, simple = FALSE))
  rendered_focal <- unlist(scene$payload$streamlines[[1]]$hits[[1]], use.names = FALSE)
  rendered_target <- unlist(scene$payload$streamlines[[1]]$targetHits[[1]], use.names = FALSE)

  expect_equal(rendered_focal, expected_focal, tolerance = 1e-6)
  expect_equal(rendered_target, expected_target, tolerance = 1e-6)
  expect_equal(scene$payload$contextLabel, "anatomical image foreground calibrated from local/QC mask")
  expect_equal(length(scene$payload$brain), 2L)
})

test_that("interactive QC viewer uses tractogram context when anatomy is absent", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "interactive-qc-viewer-context")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  fnode <- file.path(tmp, "fnode.nii.gz")
  target <- file.path(tmp, "target.nii.gz")
  context_image <- file.path(tmp, "context.nii.gz")
  tractogram <- file.path(tmp, "tracts.trk")

  fnode_arr <- array(0, dim = c(5, 5, 5))
  fnode_arr[2, 2, 2] <- 1
  target_arr <- array(0, dim = c(5, 5, 5))
  target_arr[5, 2, 2] <- 1
  context_arr <- array(0, dim = c(5, 5, 5))
  context_arr[1:5, 2, 2] <- 1
  RNifti::writeNifti(fnode_arr, fnode)
  RNifti::writeNifti(target_arr, target)
  RNifti::writeNifti(context_arr, context_image)
  ftq_write_qc_test_trk(
    tractogram,
    list(
      rbind(c(1, 2, 2), c(2, 2, 2), c(5, 2, 2)),
      rbind(c(1, 1, 1), c(1, 1, 2))
    )
  )

  manifest <- list(
    dataset_root = tmp,
    subject = "sub-test",
    description = "qc-context",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(A = basename(target)),
    coordinate = list(track_space = "voxel", discretization = "round", axis_flip = c(FALSE, FALSE, FALSE)),
    execution = list(track_reader = "stream", track_chunk_size = 2L),
    outputs = list(output_dir = file.path(tmp, "out"), write_maps = FALSE, write_report = FALSE, overwrite = TRUE)
  )
  viewer <- ftq_interactive_qc_viewer(
    manifest,
    output_html = file.path(tmp, "viewer.html"),
    sample_per_target = 1L,
    overwrite = TRUE,
    tractogram_context = context_image
  )
  scene <- ftq_scene_payload(viewer$html)
  expect_equal(scene$payload$contextLabel, "Tractogram-derived density context")
  expect_true(grepl("Tractogram-derived density context", scene$html, fixed = TRUE))
  expect_gt(length(scene$payload$brain), 0L)
})

test_that("QC point-cloud subsampling is deterministic without changing RNG state", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "interactive-qc-rng")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  background <- file.path(tmp, "background.nii.gz")
  density <- file.path(tmp, "density.nii.gz")
  RNifti::writeNifti(array(seq_len(216), dim = c(6, 6, 6)), background)
  RNifti::writeNifti(array(1, dim = c(6, 6, 6)), density)

  set.seed(812)
  seed_before <- .Random.seed
  background_cloud <- FibrTraQ:::ftq_qc_background_cloud(
    background,
    max_points = 12L
  )
  density_cloud <- FibrTraQ:::ftq_qc_density_cloud(
    density,
    max_points = 12L
  )

  expect_identical(.Random.seed, seed_before)
  expect_equal(nrow(background_cloud), 12L)
  expect_equal(nrow(density_cloud), 12L)
  expect_equal(
    background_cloud,
    FibrTraQ:::ftq_qc_background_cloud(background, max_points = 12L)
  )
  expect_equal(
    density_cloud,
    FibrTraQ:::ftq_qc_density_cloud(density, max_points = 12L)
  )
  expect_identical(.Random.seed, seed_before)
})

test_that("statistical helpers summarize, compare, and generate shifted null ROIs", {
  skip_if_not_installed("RNifti")

  result <- data.frame(
    fNode_Coords = paste(seq_len(6), 1, 1, sep = "_"),
    tracks = c(10, 12, 14, 16, 18, 20),
    local = c(1, 1, 1, 1, 1, 1),
    leave_no_target = c(1, 1, 1, 1, 1, 1),
    A = c(1, 2, 4, 8, 12, 16),
    B = c(8, 8, 7, 5, 4, 2),
    none = 0,
    stringsAsFactors = FALSE
  )
  summary <- ftq_stat_subject_summary(result, subject = "S1", group = "control", roi = "left-thalamus")
  expect_equal(summary$roi$focal_voxel_streamline_associations, sum(result$tracks))
  expect_equal(summary$targets$target[1], "A")
  expect_equal(nrow(summary$voxels), 6L)

  gradients <- ftq_stat_spatial_gradients(result, min_connected_voxels = 2L)
  expect_true(all(c("target", "axis", "slope", "q_value") %in% names(gradients)))
  expect_true(any(gradients$target == "A" & gradients$axis == "x"))

  correct <- data.frame(unit = c("u1", "u2", "u3"), assigned_tracks = c(10, 12, 14), target_hit_events = c(9, 10, 13))
  control <- data.frame(unit = c("u1", "u2", "u3"), assigned_tracks = c(10, 12, 14), target_hit_events = c(1, 2, 3))
  neg <- ftq_stat_negative_control(correct, control, by = "unit")
  expect_gt(neg$median_difference, 0)
  expect_equal(nrow(neg$paired), 3L)

  group_data <- data.frame(subject = paste0("S", 1:6), group = rep(c("A", "B"), each = 3), metric = c(3, 4, 5, 1, 2, 2))
  cmp <- ftq_stat_compare_groups(group_data, metric = "metric")
  expect_equal(cmp$metric, "metric")
  expect_true(is.finite(cmp$p_value))

  tmp <- file.path(tempdir(), "shifted-rois")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  mask <- file.path(tmp, "mask.nii.gz")
  arr <- array(0, dim = c(4, 4, 4))
  arr[2, 2, 2] <- 1
  RNifti::writeNifti(arr, mask)
  shifted <- ftq_null_shifted_rois(mask, tmp, data.frame(x = 1L, y = 0L, z = 0L), overwrite = TRUE)
  expect_true(file.exists(shifted$path))
  shifted_arr <- RNifti::readNifti(shifted$path, internal = FALSE)
  expect_equal(as.array(shifted_arr)[3, 2, 2], 1)
})
