ftq_endpoint_fixture <- function(overlap = FALSE) {
  blank <- array(0, dim = c(8, 4, 4))
  focal <- blank
  target <- blank
  context <- array(1, dim = dim(blank))
  focal[2, 2, 2] <- 1
  target[6, 2, 2] <- 1
  if (overlap) target[2, 2, 2] <- 1
  tracks <- list(
    rbind(c(2, 2, 2), c(6, 2, 2)),
    rbind(c(1, 2, 2), c(2, 2, 2), c(6, 2, 2), c(7, 2, 2)),
    rbind(c(1, 2, 2), c(2, 2, 2), c(6, 2, 2)),
    rbind(c(2, 2, 2), c(6, 2, 2), c(7, 2, 2)),
    matrix(c(2, 2, 2), nrow = 1L)
  )
  list(focal = focal, target = target, context = context, tracks = tracks)
}

ftq_endpoint_write_trk <- function(path, tracks, dims) {
  write_fixed <- function(con, text, n) {
    bytes <- rep(as.raw(0), n)
    raw <- charToRaw(text)
    if (length(raw)) bytes[seq_len(min(length(raw), n))] <- raw[seq_len(min(length(raw), n))]
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

ftq_run_endpoint_fixture <- function(backend = "r", reverse = FALSE, overlap = FALSE,
                                     endpoint_analysis = "exact_and_dilated1") {
  fixture <- ftq_endpoint_fixture(overlap = overlap)
  FibrTraQ(
    fNode = fixture$focal,
    localROI = NULL,
    cNode_list = list(target = fixture$target),
    cNode_labels = "target",
    all_tracks = fixture$tracks,
    discretization = "round",
    focal_assignment_mode = "all_focal_voxels",
    streamline_voxelization_mode = "sample_points",
    reverse_streamlines = reverse,
    track_space = "voxel",
    axis_flip = c(FALSE, FALSE, FALSE),
    backend = backend,
    endpoint_analysis = endpoint_analysis,
    endpoint_context_list = if (identical(endpoint_analysis, "none")) NULL else list(brain = fixture$context)
  )
}
