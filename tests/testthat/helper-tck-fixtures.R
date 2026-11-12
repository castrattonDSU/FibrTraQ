# Independent format writers and small spatially asymmetric fixtures. World
# points below come from explicit affine multiplication, not package mappers.
tck_fixture_write_tck <- function(path, tracks, datatype = "Float32LE") {
  endian <- if (endsWith(datatype, "BE")) "big" else "little"
  size <- if (startsWith(datatype, "Float64")) 8L else 4L
  offset <- 512L
  header <- paste0(
    "mrtrix tracks\ncount: ", length(tracks), "\n",
    "datatype: ", datatype, "\nfile: . ", offset, "\nEND\n"
  )
  bytes <- charToRaw(header)
  stopifnot(length(bytes) < offset)
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(c(bytes, raw(offset - length(bytes))), con)
  for (track in tracks) {
    writeBin(as.numeric(t(track)), con, size = size, endian = endian)
    writeBin(rep(NaN, 3L), con, size = size, endian = endian)
  }
  writeBin(rep(Inf, 3L), con, size = size, endian = endian)
  invisible(path)
}

tck_fixture_write_trk <- function(path, voxel_tracks, dims,
                                  affine_mm = diag(4L), voxel_order = "RAS") {
  spacing <- sqrt(colSums(affine_mm[1:3, 1:3, drop = FALSE]^2))
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  fixed <- function(value, n) {
    bytes <- raw(n)
    value <- charToRaw(value)
    if (length(value)) bytes[seq_along(value)] <- value
    writeBin(bytes, con)
  }
  fixed("TRACK", 6L)
  writeBin(as.integer(dims), con, size = 2L, endian = "little")
  writeBin(spacing, con, size = 4L, endian = "little")
  writeBin(rep(0, 3L), con, size = 4L, endian = "little")
  writeBin(0L, con, size = 2L, endian = "little")
  fixed("", 200L)
  writeBin(0L, con, size = 2L, endian = "little")
  fixed("", 200L)
  writeBin(as.numeric(t(affine_mm)), con, size = 4L, endian = "little")
  fixed("", 444L)
  fixed(voxel_order, 4L)
  fixed("", 4L)
  writeBin(c(1, 0, 0, 0, 1, 0), con, size = 4L, endian = "little")
  fixed("", 2L)
  writeBin(rep(0L, 6L), con, size = 1L, endian = "little")
  writeBin(as.integer(length(voxel_tracks)), con, size = 4L, endian = "little")
  writeBin(c(2L, 1000L), con, size = 4L, endian = "little")
  for (voxel in voxel_tracks) {
    # TrackVis uses physical voxel coordinates measured from a voxel corner;
    # affine_mm maps zero-based voxel centres into world millimetres.
    coords <- sweep(voxel - 0.5, 2L, spacing, "*")
    writeBin(as.integer(nrow(coords)), con, size = 4L, endian = "little")
    writeBin(as.numeric(t(coords)), con, size = 4L, endian = "little")
  }
  invisible(path)
}

tck_fixture_image <- function(values, affine_mm = diag(4L), units = "mm") {
  scale_to_mm <- switch(units, m = 1000, mm = 1, um = 0.001, 1)
  native <- affine_mm
  native[1:3, ] <- native[1:3, , drop = FALSE] / scale_to_mm
  image <- RNifti::asNifti(values, reference = list(
    sform_code = 2L, qform_code = 0L,
    pixdim = c(1, sqrt(colSums(native[1:3, 1:3, drop = FALSE]^2)), rep(1, 4L)),
    srow_x = native[1L, ], srow_y = native[2L, ], srow_z = native[3L, ]
  ))
  RNifti::pixunits(image) <- units
  image
}

tck_fixture_case <- function(affine_mm = diag(4L), units = "mm",
                             voxel_order = "RAS") {
  dims <- c(12L, 9L, 7L)
  focal <- target_a <- target_b <- array(0, dims)
  focal[2:4, 3L, 2L] <- 1
  target_a[8L, 3L, 2L] <- 1
  target_b[3L, 7L, 5L] <- 1
  voxel_tracks <- list(
    rbind(c(2, 3, 2), c(3, 3, 2), c(4, 3, 2), c(8, 3, 2)),
    rbind(c(3, 3, 2), c(3, 7, 5)),
    rbind(c(2, 3, 2), c(3, 3, 2), c(4, 3, 2)),
    rbind(c(2, 3, 2), c(1, 3, 2)),
    rbind(c(10, 7, 5), c(9, 7, 5)),
    rbind(c(1, 3, 2), c(8, 3, 2))
  )
  world_tracks <- lapply(voxel_tracks, function(voxel) {
    (cbind(voxel - 1, 1) %*% t(affine_mm))[, 1:3, drop = FALSE]
  })
  list(
    dims = dims, affine_mm = affine_mm, units = units,
    voxel_order = voxel_order, voxel_tracks = voxel_tracks,
    world_tracks = world_tracks,
    focal = tck_fixture_image(focal, affine_mm, units),
    targets = list(
      target_a = tck_fixture_image(target_a, affine_mm, units),
      target_b = tck_fixture_image(target_b, affine_mm, units)
    ),
    reference = tck_fixture_image(array(seq_len(prod(dims)), dims), affine_mm, units)
  )
}

tck_fixture_files <- function(root, fixture = tck_fixture_case(),
                              datatype = "Float32LE") {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  paths <- list(tck = file.path(root, "tracks.tck"),
                trk = file.path(root, "tracks.trk"),
                focal = file.path(root, "focal.nii.gz"),
                reference = file.path(root, "reference.nii.gz"),
                targets = stats::setNames(file.path(root,
                  paste0(names(fixture$targets), ".nii.gz")), names(fixture$targets)))
  tck_fixture_write_tck(paths$tck, fixture$world_tracks, datatype)
  tck_fixture_write_trk(paths$trk, fixture$voxel_tracks, fixture$dims,
                       fixture$affine_mm, fixture$voxel_order)
  RNifti::writeNifti(fixture$focal, paths$focal)
  RNifti::writeNifti(fixture$reference, paths$reference)
  for (name in names(fixture$targets)) {
    RNifti::writeNifti(fixture$targets[[name]], paths$targets[[name]])
  }
  paths
}

tck_fixture_expected <- function(assignment_mode = "full_streamline",
                                 voxelization = "sample_points") {
  segment <- identical(voxelization, "segment_interpolation")
  full <- identical(assignment_mode, "full_streamline")
  data.frame(
    fNode_Coords = c("2_3_2", "3_3_2", "4_3_2"),
    tracks = c(3, 3, 2) + as.integer(segment),
    local = c(1, 1, 1),
    leave_no_target = c(1, 0, 0) + as.integer(segment && !full),
    target_a = rep(1 + as.integer(segment && full), 3L),
    target_b = c(0, 1, 0), none = c(0, 0, 0),
    stringsAsFactors = FALSE
  )
}
