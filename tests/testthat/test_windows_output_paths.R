test_that("Windows output path budget includes the filename and UTF-16 terminator", {
  root <- "C:/FTQ"
  path <- paste0(root, "/", paste(rep("x", 252), collapse = ""))
  expect_equal(ftq_windows_path_length(path), 259)
  expect_invisible(ftq_check_windows_output_paths(path, root, windows = TRUE))
  expect_error(ftq_check_windows_output_paths(paste0(path, "x"), root, windows = TRUE),
    "supported limit of 259 UTF-16 code units", fixed = TRUE)
  expect_error(ftq_check_windows_output_paths(paste0(path, "x"), root, windows = TRUE),
    "output_dir (at most 5 UTF-16 code units", fixed = TRUE)
  expect_equal(ftq_windows_path_length(c("caf\u00e9", "\U0001f9e0")), c(4, 2))
  expect_error(ftq_check_windows_output_paths(paste0(substr(path, 1, 258), "\U0001f9e0"),
    root, windows = TRUE), "Longest planned path (260)", fixed = TRUE)
  expect_invisible(ftq_check_windows_output_paths(paste0(path, "longer"), root,
    windows = FALSE))
})

test_that("planned report paths retain the public names and account for options", {
  root <- "C:/FTQ/run1"
  prefix <- "S22_desc-analysis-cerebellum-regional-connectivity"
  paths <- ftq_regional_study_output_paths(root, "S22", c("cerebellum", "cerebrum"),
    "cerebellum", "analysis")
  report_dir <- file.path(root, "S22", "cerebellum", "reports", prefix)
  expect_contains(paths, file.path(report_dir, paste0(prefix, "_report.html")))
  expect_contains(paths, file.path(report_dir, "report-figures",
    paste0(prefix, "_fig_endpoint-exact-vs-dilated1.png")))
  expect_contains(paths, file.path(report_dir, "tables", paste0(prefix, "_ranked_cerebrum_voxels.csv")))
  expect_contains(paths, file.path(report_dir, "report-figures",
    paste0(prefix, "_fig_interactive-3d-qc.html")))
  short <- ftq_regional_study_output_paths(root, "S22", c("cerebellum", "cerebrum"),
    "cerebellum", "analysis", write_report = FALSE, write_maps = FALSE,
    write_exploratory_report = FALSE, write_qc_viewer = FALSE, write_endpoint = FALSE)
  expect_false(any(grepl("/reports/|/maps/|/figures/|_endpoint-", short)))
  local <- ftq_regional_study_output_paths(root, "S22", c("cerebellum", "cerebrum"),
    "cerebellum", "analysis", write_local_search = TRUE)
  expect_contains(local, file.path(report_dir, "report-figures",
    paste0(prefix, "_fig_local-search-end-voxel-distances.png")))
  long_target <- paste(rep("t", 235), collapse = "")
  labelled <- ftq_run_output_paths("C:/x", "s", long_target, write_maps = FALSE,
    write_endpoint = FALSE)
  expect_error(ftq_check_windows_output_paths(labelled, "C:/x", windows = TRUE),
    paste0("_ranked_", long_target, "_voxels.csv"), fixed = TRUE)
})

test_that("relative output roots resolve before checking without creating directories", {
  root <- tempfile("p-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  planned <- ftq_absolute_output_path(file.path("new", " # %20 caf\u00e9"))
  expect_identical(planned, file.path(normalizePath(root, winslash = "/"),
                                     "new", " # %20 caf\u00e9"))
  expect_false(dir.exists("new"))
})

test_that("file and regional workflows reject long paths before any outputs or analysis", {
  example <- ftq_example_data()
  root <- file.path(tempdir(), paste(rep("long-output-", 18), collapse = ""))
  expect_false(dir.exists(root))
  check <- ftq_check_windows_output_paths
  local_mocked_bindings(
    ftq_check_windows_output_paths = function(paths, output_dir, windows = TRUE) {
      check(paths, output_dir, windows = TRUE)
    },
    ftq_run = function(...) stop("Scientific run must not be reached"),
    .package = "FibrTraQ"
  )
  expect_error(ftq_run_files(example$tractogram, example$focal_mask,
    example$target_masks, root, subject = "subject", verbose = FALSE),
    "Windows output path exceeds", fixed = TRUE)
  expect_false(dir.exists(root))
  expect_error(ftq_run_regional_study(c(subject = example$tractogram),
    list(subject = c(focal = example$focal_mask, example$target_masks)),
    root, focal_regions = "focal", verbose = FALSE),
    "Windows output path exceeds", fixed = TRUE)
  expect_false(dir.exists(root))
})

test_that("manifest workflow rejects long paths before scientific preflight", {
  example <- ftq_example_data()
  root <- file.path(tempdir(), paste(rep("long-output-", 18), collapse = ""))
  check <- ftq_check_windows_output_paths
  local_mocked_bindings(
    ftq_check_windows_output_paths = function(paths, output_dir, windows = TRUE) {
      check(paths, output_dir, windows = TRUE)
    },
    ftq_preflight = function(...) stop("Scientific preflight must not be reached"),
    .package = "FibrTraQ"
  )
  manifest <- list(subject = "subject", description = "test",
    tractogram = example$tractogram, fnode = example$focal_mask,
    targets = as.list(example$target_masks), outputs = list(output_dir = root))
  expect_error(ftq_run(manifest), "Windows output path exceeds", fixed = TRUE)
  expect_false(dir.exists(root))
  short_root <- tempfile("p-")
  manifest$outputs$output_dir <- short_root
  manifest$outputs$write_report <- FALSE
  manifest$report <- list(context_source = "tractogram", tractogram_context = list(
    output_dir = root, prefix = "density"))
  expect_error(ftq_run(manifest), "tractogram-density.nii.gz", fixed = TRUE)
  expect_false(dir.exists(short_root))
  expect_false(dir.exists(root))
})
