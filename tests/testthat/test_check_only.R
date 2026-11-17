test_that("check-only checks both formats and every focal without writing results", {
  root <- tempfile("fc")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- tck_fixture_files(file.path(root, "in"))
  masks <- c(A = files$focal, B = unname(files$targets[1L]))
  output <- file.path(root, "out")
  args <- list(tractograms = c(one = files$trk, two = files$tck),
               region_masks = list(one = masks, two = masks), local_masks = "union",
               output_dir = output, check_only = TRUE, verbose = FALSE, open = TRUE)
  testthat::local_mocked_bindings(
    ftq_run = function(...) stop("Unexpected counting"),
    ftq_qc_open_output = function(...) stop("Unexpected viewer"),
    ftq_write_report_page = function(...) stop("Unexpected report"),
    .package = "FibrTraQ"
  )
  check <- do.call(fibrtraq_study, args)
  expect_true(check$passed)
  expect_identical(check$index$subject, rep(c("one", "two"), each = 2))
  expect_identical(check$index$focal_region, rep(c("A", "B"), 2))
  expect_true(all(check$index$passed))
  expect_true(any(check$checks$check == "tck_world_coordinate_mapping"))
  expect_true(any(check$checks$check == "local_readable"))
  expect_false(file.exists(output))
  expect_output(print(check), "Review anatomical alignment independently")

  # A blocked output path is local to this subject and checking preserves the file.
  dir.create(output)
  writeLines("keep", file.path(output, "one"))
  args$overwrite <- TRUE
  check <- do.call(fibrtraq_study, args)
  expect_false(check$passed)
  expect_identical(check$index$passed, c(FALSE, FALSE, TRUE, TRUE))
  expect_true(any(check$checks$check == "output_directory"))
  expect_identical(readLines(file.path(output, "one")), "keep")
  expect_identical(list.files(output), "one")
})

test_that("input checks identify independent header and grid errors by subject", {
  root <- tempfile("fc")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- tck_fixture_files(file.path(root, "in"))
  malformed <- file.path(root, "bad.trk")
  writeLines("invalid", malformed)
  shifted <- file.path(root, "off.nii.gz")
  affine <- diag(4)
  affine[1, 4] <- 4
  image <- RNifti::readNifti(files$focal)
  RNifti::writeNifti(tck_fixture_image(image, affine), shifted)
  masks <- c(A = files$focal, B = unname(files$targets[1L]))
  off_grid <- c(A = files$focal, B = shifted)
  output <- file.path(root, "out")
  check <- fibrtraq_study(
    tractograms = c(good = files$trk, header = malformed, grid = files$tck,
                    mapping = files$tck),
    region_masks = list(good = masks, header = masks, grid = off_grid, mapping = masks),
    axis_flip = list(good = NULL, header = NULL, grid = NULL,
                     mapping = c(TRUE, FALSE, FALSE)),
    focal_regions = "A", output_dir = output, check_only = TRUE, verbose = FALSE
  )
  expect_identical(check$index$passed, c(TRUE, FALSE, FALSE, FALSE))
  errors <- subset(check$checks, severity == "error" & !passed)
  expect_true(any(errors$subject == "header" & errors$check == "trk_header_readable"))
  expect_true(any(errors$subject == "grid" & grepl("off.nii.gz", errors$detail, fixed = TRUE)))
  expect_true(any(errors$subject == "mapping" & errors$check == "tck_world_coordinate_mapping"))
  expect_false(file.exists(output))
})

test_that("single and equivalent one-subject studies use identical defaults and counts", {
  root <- tempfile("fc")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- tck_fixture_files(file.path(root, "in"))
  masks <- c(focal = files$focal, files$targets)
  settings <- list(endpoint_analysis = "none", write_maps = FALSE, write_report = FALSE,
                   write_exploratory_report = FALSE, write_qc_viewer = FALSE,
                   verbose = FALSE)
  one <- c(list(tractogram = files$tck, focal_mask = files$focal, target_masks = files$targets,
                subject = "001", focal_label = "focal", output_dir = file.path(root, "s")), settings)
  cohort <- c(list(tractograms = c("001" = files$tck), region_masks = list("001" = masks),
                   focal_regions = "focal", output_dir = file.path(root, "c")), settings)
  single_check <- do.call(fibrtraq_single, c(one, list(check_only = TRUE)))
  cohort_check <- do.call(fibrtraq_study, c(cohort, list(check_only = TRUE)))
  expect_true(single_check$passed)
  expect_identical(cohort_check$index, single_check$index)
  expect_false(file.exists(one$output_dir))
  expect_false(file.exists(cohort$output_dir))

  one$verbose <- TRUE
  messages <- capture.output(single <- do.call(fibrtraq_single, one), type = "message")
  expect_true(any(grepl(single$report_index, messages, fixed = TRUE)))
  expect_true(any(grepl(paste0("Share this folder: ", single$output_dir), messages, fixed = TRUE)))
  expect_message(study <- do.call(fibrtraq_study, cohort), NA)
  for (field in c("result_csv", "rsc_csv", "fate_csv")) {
    expect_identical(readLines(study$index[[field]]), readLines(single$index[[field]]))
  }
  expect_false(dir.exists(file.path(study$output_dir, "001", "prepared")))
  expect_output(print(study), "Open this report:")

  # A rerun check never overwrites an existing report, including with overwrite=TRUE.
  hash <- tools::md5sum(study$report_index)
  expect_false(do.call(fibrtraq_study, c(cohort, list(check_only = TRUE)))$passed)
  expect_true(do.call(fibrtraq_study, c(cohort, list(check_only = TRUE, overwrite = TRUE)))$passed)
  expect_identical(tools::md5sum(study$report_index), hash)
})

test_that("case-colliding target names cannot overwrite files", {
  root <- tempfile("fc")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- tck_fixture_files(file.path(root, "in"))
  output <- file.path(root, "out")
  dir.create(output)
  previous <- file.path(output, "fibrtraq-run-index.csv")
  writeLines("prior outputs", previous)
  unchanged <- tools::md5sum(c(files$targets, files$focal, previous))
  targets <- stats::setNames(files$targets, c("Pons", "pons"))
  one <- list(tractogram = files$tck, focal_mask = files$focal, target_masks = targets,
              focal_label = "A", subject = "s")
  cohort <- list(tractograms = c(s = files$tck), region_masks = list(s = c(A = files$focal, targets)),
                 focal_regions = "A")
  for (check_only in c(TRUE, FALSE)) {
    settings <- list(check_only = check_only, overwrite = TRUE, output_dir = output, verbose = FALSE)
    expect_error(do.call(fibrtraq_single, c(one, settings)), "conflicting regions: Pons, pons", fixed = TRUE)
    expect_error(do.call(fibrtraq_study, c(cohort, settings)), "conflicting regions: Pons, pons", fixed = TRUE)
  }
  preflight <- ftq_preflight(list(dataset_root = root, subject = "s", description = "test", fnode = files$focal,
                                  tractogram = files$tck, targets = as.list(targets),
                                  outputs = list(output_dir = output, overwrite = TRUE)))
  duplicates <- subset(preflight$checks, check == "target_map_names_unique")
  expect_identical(duplicates$passed, FALSE)
  expect_match(duplicates$detail, "Pons, pons", fixed = TRUE)
  expect_identical(tools::md5sum(names(unchanged)), unchanged)
  expect_identical(list.files(output), basename(previous))
})

test_that("existing files for requested report stages are detected before counting", {
  root <- tempfile("fc")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- tck_fixture_files(file.path(root, "in"))
  output <- file.path(root, "out")
  focal_dir <- file.path(output, "001", "A")
  figures <- file.path(focal_dir, "figures")
  dir.create(figures, recursive = TRUE)
  prefix <- ftq_prefix("001", desc = "analysis-A-regional-connectivity")
  orphan <- file.path(figures, paste0(prefix, "_fig_fate-summary.png"))
  writeLines("prior figure", orphan)
  args <- list(tractogram = files$tck, focal_mask = files$focal, target_masks = files$targets,
               subject = "001", focal_label = "A", output_dir = output, verbose = FALSE)
  run <- function(check_only = FALSE) do.call(fibrtraq_single, c(args, list(check_only = check_only)))
  check <- run(TRUE)
  expect_false(check$passed)
  expect_true(any(grepl(basename(orphan), check$checks$detail, fixed = TRUE)))
  expect_error(run(), "already exist", fixed = TRUE)
  expect_identical(readLines(orphan), "prior figure")
  expect_false(file.exists(file.path(output, "fibrtraq-run-index.csv")))

  args$write_report <- FALSE
  expect_true(run(TRUE)$passed)
  args$write_report <- TRUE
  args$overwrite <- TRUE
  expect_true(run(TRUE)$passed)
  # Another requested output's filename is occupied by a directory.
  occupied <- file.path(figures, paste0(prefix, "_fig_rsc-summary.png"))
  dir.create(occupied)
  expect_false(run(TRUE)$passed)
  expect_error(run(), "cannot be replaced", fixed = TRUE)
  expect_true(dir.exists(occupied))
  expect_false(file.exists(file.path(output, "fibrtraq-run-index.csv")))
})

test_that("a target cannot replace the built-in ranked table through a case variant", {
  root <- tempfile("fc")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- tck_fixture_files(file.path(root, "in"))
  output <- file.path(root, "out")
  table_dir <- file.path(output, "s", "focal", "tables")
  dir.create(table_dir, recursive = TRUE)
  existing <- ftq_ranked_table_path(table_dir, "s_desc-analysis-focal-regional-connectivity", "leave_no_target")
  writeLines("prior ranks", existing)
  targets <- c(LEAVE_NO_TARGET = unname(files$targets[1L]))
  one <- list(tractogram = files$tck, focal_mask = files$focal, target_masks = targets,
              subject = "s", focal_label = "focal")
  cohort <- list(tractograms = c(s = files$tck), region_masks = list(s = c(focal = files$focal, targets)))
  for (check_only in c(TRUE, FALSE)) {
    settings <- list(check_only = check_only, output_dir = output, verbose = FALSE, overwrite = TRUE)
    expect_error(do.call(fibrtraq_single, c(one, settings)), "ranked tables: LEAVE_NO_TARGET", fixed = TRUE)
    expect_error(do.call(fibrtraq_study, c(cohort, settings)), "ranked tables: LEAVE_NO_TARGET", fixed = TRUE)
  }
  # A region used only as focal produces no target-ranked table of this name.
  expect_true(do.call(fibrtraq_study, c(cohort, list(check_only = TRUE, output_dir = output,
                      focal_regions = "LEAVE_NO_TARGET", verbose = FALSE, overwrite = TRUE)))$passed)
  internal <- ftq_preflight(list(dataset_root = root, subject = "s", description = "test",
                                 fnode = files$focal, tractogram = files$tck, targets = as.list(targets),
                                 outputs = list(output_dir = output, overwrite = TRUE)))
  expect_identical(subset(internal$checks, check == "target_labels_not_reserved")$passed, FALSE)
  expect_identical(readLines(existing), "prior ranks")
  expect_length(list.files(output, recursive = TRUE), 1L)
})
