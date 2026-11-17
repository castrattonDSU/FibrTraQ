test_that("installed validation manifest templates parse", {
  manifest_dir <- file.path(archive_root, "validation", "manifests")
  skip_if(!nzchar(manifest_dir), "validation templates not installed")

  manifests <- list.files(manifest_dir, pattern = "\\.yml$", full.names = TRUE)
  expect_gt(length(manifests), 0)
  parsed <- lapply(manifests, ftq_read_manifest)

  expect_true(all(vapply(parsed, function(x) length(x$target_labels) > 0L, logical(1))))
  expect_true(all(vapply(parsed, function(x) identical(x$track_space, "auto"), logical(1))))
  expect_true(all(vapply(parsed, function(x) isTRUE(x$write_report), logical(1))))
})

test_that("validation suite summarizes manifest preflight rows", {
  skip_if_not_installed("RNifti")

  tmp <- tempdir()
  fnode <- file.path(tmp, "suite-fnode.nii.gz")
  target <- file.path(tmp, "suite-target.nii.gz")
  tractogram <- file.path(tmp, "suite-tracts.trk")
  RNifti::writeNifti(array(c(1, rep(0, 26)), dim = c(3, 3, 3)), fnode)
  RNifti::writeNifti(array(c(0, 1, rep(0, 25)), dim = c(3, 3, 3)), target)
  file.create(tractogram)

  manifest <- list(
    dataset_root = tmp,
    subject = "sub-suite",
    description = "validation-suite",
    fnode = basename(fnode),
    tractogram = basename(tractogram),
    targets = list(T = basename(target)),
    outputs = list(output_dir = file.path(tmp, "suite-out"), overwrite = FALSE)
  )
  manifest_path <- file.path(tmp, "suite-manifest.yml")
  yaml::write_yaml(manifest, manifest_path)
  summary_csv <- file.path(tmp, "suite-summary.csv")

  summary <- ftq_run_validation_suite(
    manifests = manifest_path,
    output_csv = summary_csv,
    preflight_only = TRUE
  )

  expect_true(file.exists(summary_csv))
  expect_equal(nrow(summary), 1L)
  expect_equal(summary$run_status, "preflight_only")
  expect_true(summary$preflight_passed)
  expect_equal(summary$preflight_recommendations, 0L)
  expect_true(is.na(summary$validation_passed))
})

test_that("ISMRM 2015 manifest generator maps configured endpoints", {
  tmp <- file.path(tempdir(), "ismrm-generator")
  scoring <- file.path(tmp, "scoring")
  manifest_dir <- file.path(tmp, "manifests")
  dir.create(file.path(scoring, "bundles"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(scoring, "bundles", "sub_bundles"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(scoring, "ROI", "endpoints"), recursive = TRUE, showWarnings = FALSE)

  file.create(file.path(scoring, "bundles", "CP.trk"))
  file.create(file.path(scoring, "bundles", "sub_bundles", "CC_u_shaped.trk"))
  file.create(file.path(scoring, "ROI", "endpoints", "CP_head.nii.gz"))
  file.create(file.path(scoring, "ROI", "endpoints", "CP_tail.nii.gz"))
  file.create(file.path(scoring, "ROI", "endpoints", "CC_left.nii.gz"))
  file.create(file.path(scoring, "ROI", "endpoints", "CC_right.nii.gz"))
  file.create(file.path(scoring, "t1.nii.gz"))
  jsonlite::write_json(
    list(
      CP = list(
        all_mask = "ROI/all_masks/CP.nii.gz",
        head = "ROI/endpoints/CP_head.nii.gz",
        tail = "ROI/endpoints/CP_tail.nii.gz"
      ),
      CC_u_shaped = list(
        all_mask = "ROI/all_masks/CC_u_shaped.nii.gz",
        head = "ROI/endpoints/CC_left.nii.gz",
        tail = "ROI/endpoints/CC_right.nii.gz"
      )
    ),
    file.path(scoring, "config_file_segmentation.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  generated <- ftq_make_ismrm2015_manifests(
    scoring_dir = scoring,
    output_dir = manifest_dir,
    run_output_root = file.path(tmp, "outputs"),
    force = TRUE
  )

  expect_equal(generated$bundle, c("CP", "CC_u_shaped"))
  expect_true(all(file.exists(generated$manifest_path)))

  parsed <- ftq_read_manifest(generated$manifest_path[1])
  expect_equal(parsed$description, "ismrm2015-CP-groundtruth")
  expect_true(grepl("bundles/CP.trk$", ftq_test_path(parsed$tractogram)))
  expect_true(grepl("ROI/endpoints/CP_head.nii.gz$", ftq_test_path(parsed$fnode)))
  expect_equal(parsed$target_labels, "tail")
  expect_true(grepl("ROI/endpoints/CP_tail.nii.gz$", ftq_test_path(parsed$target_paths[1])))
  expect_true(grepl("t1.nii.gz$", ftq_test_path(parsed$background_image)))
})

test_that("ISMRM 2015 control manifest generators build specificity tests", {
  skip_if_not_installed("RNifti")

  tmp <- file.path(tempdir(), "ismrm-control-generator")
  scoring <- file.path(tmp, "scoring")
  dir.create(file.path(scoring, "bundles"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(scoring, "ROI", "endpoints"), recursive = TRUE, showWarnings = FALSE)

  file.create(file.path(scoring, "bundles", "CP.trk"))
  file.create(file.path(scoring, "bundles", "MCP.trk"))
  RNifti::writeNifti(array(c(1, rep(0, 26)), dim = c(3, 3, 3)),
                     file.path(scoring, "ROI", "endpoints", "CP_head.nii.gz"))
  RNifti::writeNifti(array(c(0, 1, rep(0, 25)), dim = c(3, 3, 3)),
                     file.path(scoring, "ROI", "endpoints", "CP_tail.nii.gz"))
  RNifti::writeNifti(array(c(0, 0, 1, rep(0, 24)), dim = c(3, 3, 3)),
                     file.path(scoring, "ROI", "endpoints", "MCP_head.nii.gz"))
  RNifti::writeNifti(array(c(0, 0, 0, 1, rep(0, 23)), dim = c(3, 3, 3)),
                     file.path(scoring, "ROI", "endpoints", "MCP_tail.nii.gz"))
  RNifti::writeNifti(array(seq_len(27), dim = c(3, 3, 3)),
                     file.path(scoring, "t1.nii.gz"))
  jsonlite::write_json(
    list(
      CP = list(head = "ROI/endpoints/CP_head.nii.gz", tail = "ROI/endpoints/CP_tail.nii.gz"),
      MCP = list(head = "ROI/endpoints/MCP_head.nii.gz", tail = "ROI/endpoints/MCP_tail.nii.gz")
    ),
    file.path(scoring, "config_file_segmentation.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  wrong <- ftq_make_ismrm2015_wrong_target_manifests(
    scoring_dir = scoring,
    output_dir = file.path(tmp, "wrong"),
    bundles = "CP",
    target_bundles = c("CP", "MCP"),
    run_output_root = file.path(tmp, "wrong-out"),
    force = TRUE
  )
  expect_equal(nrow(wrong), 1L)
  expect_equal(wrong$target_bundle, "MCP")
  wrong_parsed <- ftq_read_manifest(wrong$manifest_path)
  expect_equal(wrong_parsed$target_labels, "wrong_MCP")
  expect_false(wrong_parsed$write_maps)

  wrong_matrix <- ftq_make_ismrm2015_wrong_target_manifests(
    scoring_dir = scoring,
    output_dir = file.path(tmp, "wrong-matrix"),
    bundles = "CP",
    target_bundles = c("CP", "MCP"),
    run_output_root = file.path(tmp, "wrong-matrix-out"),
    combine_targets = TRUE,
    force = TRUE
  )
  expect_equal(nrow(wrong_matrix), 1L)
  expect_equal(wrong_matrix$target_count, 1L)
  wrong_matrix_parsed <- ftq_read_manifest(wrong_matrix$manifest_path)
  expect_equal(wrong_matrix_parsed$target_labels, "wrong_MCP")
  expect_true(grepl("all-wrong-targets$", wrong_matrix_parsed$output_dir))

  shifted <- ftq_make_ismrm2015_shifted_target_manifests(
    scoring_dir = scoring,
    output_dir = file.path(tmp, "shifted"),
    shifted_roi_dir = file.path(tmp, "shifted-rois"),
    shift = c(1L, 0L, 0L),
    bundles = "CP",
    run_output_root = file.path(tmp, "shifted-out"),
    force = TRUE
  )
  expect_equal(nrow(shifted), 1L)
  expect_true(file.exists(shifted$target))
  shifted_arr <- RNifti::readNifti(shifted$target, internal = FALSE)
  expect_equal(sum(as.array(shifted_arr) != 0), 1)
  shifted_parsed <- ftq_read_manifest(shifted$manifest_path)
  expect_equal(shifted_parsed$target_labels, "shifted_tail")
  expect_true(grepl("shifted-rois", shifted_parsed$target_paths[1]))

  submission <- file.path(tmp, "submission.trk")
  file.create(submission)
  sweep <- ftq_make_ismrm2015_submission_sweep_manifests(
    scoring_dir = scoring,
    tractograms = submission,
    output_dir = file.path(tmp, "sweep"),
    bundles = c("CP", "MCP"),
    run_output_root = file.path(tmp, "sweep-out"),
    force = TRUE
  )
  expect_equal(nrow(sweep), 2L)
  expect_equal(sweep$bundle, c("CP", "MCP"))
  sweep_parsed <- ftq_read_manifest(sweep$manifest_path[1])
  expect_true(grepl("submission.trk$", sweep_parsed$tractogram))
})
