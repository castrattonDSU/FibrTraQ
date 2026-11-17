test_that("a target NIfTI counts each qualifying streamline once per focal voxel", {
  root <- tempfile("tm-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dims <- c(10L, 8L, 5L)
  a <- c(2, 3, 2)
  b <- c(3, 3, 2)
  c0 <- c(4, 3, 2)
  untouched <- c(5, 6, 4)
  m1 <- c(7, 3, 2)
  m2 <- c(8, 3, 2)
  p <- c(2, 6, 2)
  focal <- medullaspinal <- pontine <- array(0, dims)
  focal[rbind(a, b, c0, untouched)] <- 1
  medullaspinal[rbind(m1, m2)] <- 1
  pontine[rbind(p)] <- 1
  tracks <- list(
    rbind(a, b, a, b, m1, m2, m1),          # revisits two focal and two target voxels
    rbind(a, p),                             # pontine only
    rbind(b, m2, c(8, 6, 2), p, c(3, 6, 2)), # both targets, neither is an endpoint
    rbind(a, b, b, a),                       # focal only
    rbind(m2, p),                            # target only at stored points
    rbind(c(1, 3, 2), c(9, 3, 2)),           # crosses focal and target between samples
    rbind(c0, c(4, 5, 2))                    # focal only
  )
  fixtures <- list(dims = dims, affine_mm = diag(4), voxel_order = "RAS",
                   voxel_tracks = tracks, world_tracks = lapply(tracks, function(v) v - 1),
                   focal = tck_fixture_image(focal),
                   targets = list(medullaspinal = tck_fixture_image(medullaspinal),
                                  pontine = tck_fixture_image(pontine)),
                   reference = tck_fixture_image(array(1, dims)))
  paths <- tck_fixture_files(file.path(root, "in"), fixtures)

  check <- function(format, backend, interpolate = FALSE) {
    output <- file.path(root, paste0(format, backend, as.integer(interpolate)))
    check_report <- format == "tck" && backend == "cpp"
    args <- list(tractogram = paths[[format]], focal_mask = paths$focal,
                 target_masks = paths$targets, reference_image = paths$reference,
                 output_dir = output, subject = "s", focal_label = "f",
                 track_space = if (format == "trk") "trackvis" else "world",
                 axis_flip = rep(FALSE, 3L), track_chunk_size = 2L, backend = backend,
                 write_report = check_report, write_exploratory_report = FALSE, write_qc_viewer = FALSE,
                 verbose = FALSE, progress = FALSE)
    all_local <- format == "trk" && backend == "cpp"
    if (all_local) args$local_mask <- paths$reference  # the reference is all ones
    if (interpolate) args$streamline_voxelization_mode <- "segment_interpolation"
    run <- do.call(fibrtraq_single, args)
    expect_true(run$index$validation_passed)
    if (check_report) {
      html <- paste(readLines(run$report_index, warn = FALSE), collapse = "\n")
      sampling <- if (interpolate) "interpolated" else "stored"
      expect_match(html, paste("each assigned focal voxel from its", sampling, "streamline points"),
                   fixed = TRUE)
    }
    if (all_local) {
      counts <- read.csv(run$index$result_csv)
      expect_identical(counts$local, counts$tracks)
    }
    # Whole-volume oracle (including zeros outside the focal mask) is assigned
    # from the synthetic paths above, independently of the mapper and CSV.
    m <- q <- array(0, dims)
    m[rbind(a, b, c0)] <- if (interpolate) c(2, 3, 3) else c(1, 2, 0)
    q[rbind(a, b, c0)] <- if (interpolate) c(1, 1, 1) else c(1, 1, 0)
    for (label in names(paths$targets)) {
      map <- list.files(file.path(output, "s", "f", "maps"),
                        paste0("_target-", label, "[.]nii[.]gz$"), full.names = TRUE)
      expect_length(map, 1L)
      image <- as.array(RNifti::readNifti(map))
      expected <- if (label == "medullaspinal") m else q
      expect_identical(dim(image), dim(expected))
      expect_identical(as.numeric(image), as.numeric(expected))
    }
  }
  for (format in c("trk", "tck")) {
    check(format, "r")
    check(format, "cpp")
  }
  check("tck", "cpp", interpolate = TRUE)
})

test_that("counts above float32 integer precision are saved exactly", {
  root <- tempfile("tm-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  template <- file.path(root, "f.nii.gz")
  RNifti::writeNifti(tck_fixture_image(array(1, c(2, 2, 2))), template)
  counts <- c(2^24 + 1, 2^24 + 3)
  result <- data.frame(fNode_Coords = c("1_1_1", "2_1_1"), tracks = counts,
                       local = 0, leave_no_target = 0, medullaspinal = counts)
  paths <- FibrTraQ:::ftq_write_maps(result, template, root, "s", "medullaspinal")
  map <- RNifti::readNifti(paths[["target-medullaspinal"]])
  expect_identical(as.numeric(map), c(counts, rep(0, 6)))
  expect_equal(RNifti::niftiHeader(paths[["tracks-binary"]])$datatype, 16L)
})
