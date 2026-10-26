public_exploratory_fixture <- function() {
  data.frame(
    fNode_Coords = c("1_1_1", "2_2_2", "3_3_3", "1_2_1", "2_1_2"),
    tracks = c(4L, 2L, 0L, 1L, 1L),
    local = c(1L, 0L, 0L, 1L, 0L),
    leave_no_target = c(0L, 0L, 0L, 0L, 1L),
    A = c(2L, 2L, 0L, 0L, 0L),
    B = c(2L, 0L, 0L, 0L, 0L), none = c(0L, 0L, 1L, 0L, 0L)
  )
}

test_that("CSV exploratory reports recover matching relocated method metadata", {
  root <- tempfile("exploratory-provenance-")
  dir.create(root)
  input <- file.path(root, "S1c_fibrtraq_per-voxel-hits.csv")
  meta <- file.path(root, "S1c_fibrtraq_run-metadata.json")
  result <- public_exploratory_fixture()
  utils::write.csv(result, input, row.names = FALSE)
  for (mode in c("full_streamline", "endpoint_voxels", "local_search")) {
    # The metadata's original workstation path may differ after relocation.
    jsonlite::write_json(list(
      files = list(result_csv = "C:\\old folder\\S1c_fibrtraq_per-voxel-hits.csv"),
      manifest = list(subject = "S1c", target_intersection_mode = mode,
        local_search_radius_mm = 3.25, focal_assignment_mode = "first_intersection"),
      summary = list(target_intersection_mode = mode, local_search_radius_mm = 3.25)
    ), meta, auto_unbox = TRUE)
    from_csv <- ftq_exploratory_summary(input, write_maps = FALSE)
    direct <- ftq_exploratory_summary(result, assignment_mode = mode,
      search_radius_mm = 3.25, write_maps = FALSE)
    expect_equal(from_csv$data, direct$data)
    expect_equal(from_csv$target_totals, direct$target_totals)
    expect_identical(from_csv$provenance$target_intersection_mode, mode)
    expect_equal(from_csv$provenance$local_search_radius_mm,
      if (mode == "local_search") 3.25 else NA_real_)
    context <- attr(from_csv$provenance, "report_context")
    expect_identical(context$mode_source, "Companion run metadata")
    expect_identical(context$subject, "S1c")
  }
  expect_error(ftq_exploratory_summary(input, assignment_mode = "full_streamline"),
    "conflicts with recorded")
  expect_error(ftq_exploratory_summary(input, search_radius_mm = 2),
    "conflicts with recorded")
  jsonlite::write_json(list(
    summary = list(target_intersection_mode = "local_search"),
    manifest = list(target_intersection_mode = "full_streamline")
  ), meta, auto_unbox = TRUE)
  expect_error(ftq_exploratory_summary(input), "Conflicting target-assignment modes")
  jsonlite::write_json(list(files = list(result_csv = "other.csv")), meta, auto_unbox = TRUE)
  expect_error(ftq_exploratory_summary(input), "different result CSV")
})

test_that("missing CSV method provenance is explicit while preserving legacy values", {
  input <- tempfile(fileext = ".csv")
  result <- public_exploratory_fixture()
  utils::write.csv(result, input, row.names = FALSE)
  expect_warning(default <- ftq_exploratory_summary(input, write_maps = FALSE),
    "not a verified method")
  explicit <- ftq_exploratory_summary(input, assignment_mode = "full_streamline", write_maps = FALSE)
  expect_equal(default$data, explicit$data)
  expect_equal(default$target_totals, explicit$target_totals)
  expect_identical(default$provenance$target_intersection_mode, "full_streamline")
  expect_match(attr(default$provenance, "report_context")$mode_source, "unverified")
  expect_warning(local <- ftq_exploratory_summary(input, assignment_mode = "local_search"),
    "2 mm default is unverified")
  expect_identical(local$provenance$local_search_radius_mm, 2)
  expect_equal(explicit$data$primary_target, c("A", "A", "none", "none", "none"))
  expect_equal(explicit$data$primary_target_fraction, c(0.5, 1, 0, 0, 0))
  expect_equal(explicit$data$target_entropy_fraction, c(1, 0, 0, 0, 0))
  expect_equal(explicit$data$leave_no_target_fraction, c(0, 0, 0, 0, 1))
  expect_identical(names(explicit$provenance), c("target_intersection_mode", "local_search_radius_mm"))
  root <- tempfile("unverified-exploratory-")
  dir.create(root)
  output <- file.path(root, "report.html")
  FibrTraQ:::ftq_write_exploratory_html_report(output, default, list())
  html <- paste(readLines(output), collapse = "")
  expect_match(html, "Method provenance requires attention", fixed = TRUE)
  expect_match(html, "not a verified method for these data", fixed = TRUE)
  expect_false(grepl("Output definitions|Research maps", html))
  expect_identical(list.files(root), "report.html")
})

test_that("target report writes only seven public artifacts with portable links", {
  skip_if_not_installed("RNifti")
  root <- tempfile("explore # 100%20 ")
  dir.create(root)
  template <- file.path(root, "template.nii.gz")
  image <- RNifti::asNifti(array(0, c(3, 3, 3)))
  RNifti::pixdim(image) <- c(1.5, 2, 3)
  RNifti::writeNifti(image, template)
  FibrTraQ:::ftq_write_output_guide(root)
  result <- public_exploratory_fixture()
  attr(result, "fibrtraq_target_intersection_mode") <- "full_streamline"
  attr(result, "fibrtraq_focal_assignment_mode") <- "first_intersection"
  output_dir <- file.path(root, "report # 100%20")
  plotted <- list()
  montage <- FibrTraQ:::ftq_plot_montage
  categorical <- FibrTraQ:::ftq_plot_categorical_montage
  local_mocked_bindings(
    ftq_plot_montage = function(arr, path, ...) {
      plotted[[basename(path)]] <<- arr
      montage(arr, path, ...)
    },
    ftq_plot_categorical_montage = function(arr, lut, path, ...) {
      plotted[[basename(path)]] <<- arr
      categorical(arr, lut, path, ...)
    }, .package = "FibrTraQ"
  )
  report <- ftq_exploratory_report(result, output_dir,
    template, prefix = "S1c' test", report_slices = c(1, 3))
  html <- paste(readLines(report$html, warn = FALSE), collapse = "\n")
  expect_match(html, "Median primary-target fraction.*all focal voxels</td><td>0</td>")
  expect_match(html, "Median primary-target fraction.*target-positive voxels</td><td>0.75</td>")
  expect_match(html, "Voxels with target associations</td><td>2 / 5 \\(40%\\)</td>")
  expect_match(html, "Zero-assignment focal voxels</td><td>1 / 5 \\(20%\\)</td>")
  expect_match(html, "Equal maxima select the first target", fixed = TRUE)
  expect_match(html, "configured local ROI", fixed = TRUE)
  expect_match(html, "first_intersection", fixed = TRUE)
  expect_match(html, "FibrTraQ target report", fixed = TRUE)
  expect_false(grepl("Advanced maps|Research maps|exploratory-per-voxel.csv|exploratory-provenance.csv", html))
  expect_match(html, "<details><summary>Target association table", fixed = TRUE)
  expect_match(html, "Output definitions and interpretation guide", fixed = TRUE)
  expect_false(grepl(" onerror='", html, fixed = TRUE))
  expect_match(html, "%27", fixed = TRUE)
  hrefs <- regmatches(html, gregexpr("(?:href|src)='[^']+'", html, perl = TRUE))[[1L]]
  paths <- sub("^[^=]+='(.*)'$", "\\1", hrefs)
  paths <- paths[!grepl("^#", paths)]
  expect_true(all(file.exists(file.path(dirname(report$html), utils::URLdecode(paths)))))
  expect_false(any(grepl("%5C|^file:", paths)))
  # Compare existing filesystem paths independently of Windows separators.
  canonical <- function(paths) normalizePath(paths, winslash = "/", mustWork = TRUE)
  expect_identical(canonical(dirname(report$html)), canonical(file.path(output_dir, "figures")))
  expect_identical(basename(report$html), "S1c' test_fibrtraq_target-report.html")
  expect_identical(canonical(dirname(report$summary$outputs$lut_csv)), canonical(output_dir))
  expect_identical(sort(canonical(list.files(output_dir, recursive = TRUE, full.names = TRUE))),
    sort(canonical(c(report$html, unlist(report$figures, use.names = FALSE), report$summary$outputs$lut_csv))))
  expect_length(list.files(output_dir, recursive = TRUE), 7L)
  expect_length(report$summary$outputs$maps, 0L)
  expect_null(report$summary$outputs$summary_csv)
  expect_null(report$summary$outputs$target_totals_csv)
  expect_null(report$summary$outputs$provenance_csv)
  expect_error(ftq_exploratory_report(result, output_dir, template,
    prefix = "S1c' test", overwrite = FALSE), "already exist")

  # Explicit research exports retain all original metrics and geometry. The
  # public figures receive the same values without any intermediate map files.
  research <- ftq_exploratory_summary(result, file.path(root, "research"),
    template_image = template, prefix = "research")
  expect_length(research$outputs$maps, 10L)
  expect_true(all(file.exists(unlist(research$outputs[1:4]))))
  csv <- utils::read.csv(research$outputs$summary_csv, check.names = FALSE)
  expect_equal(csv, report$summary$data)
  expect_equal(research$target_totals, report$summary$target_totals)
  expect_identical(names(csv), c("fNode_Coords", "x", "y", "z", "tracks", "local",
    "leave_no_target", "leave_with_target", "total_target_hits", "target_hits_per_assigned_track",
    "target_diversity", "primary_target", "primary_target_index", "primary_target_hits",
    "primary_target_fraction", "target_entropy_bits", "target_entropy_fraction",
    "leave_no_target_fraction", "leave_with_target_fraction"))
  plotted_names <- c(primary_target_index = basename(report$figures$primary_target),
    target_entropy_fraction = basename(report$figures$entropy),
    target_diversity = basename(report$figures$diversity),
    leave_no_target_fraction = basename(report$figures$no_target))
  for (name in names(research$outputs$maps)) {
    arr <- RNifti::readNifti(research$outputs$maps[[name]], internal = FALSE)
    expect_equal(dim(arr), c(3L, 3L, 3L))
    expect_equal(RNifti::xform(arr), RNifti::xform(image), tolerance = 1e-7)
    expect_equal(RNifti::pixdim(arr), RNifti::pixdim(image))
    expect_equal(as.numeric(arr[as.matrix(csv[, c("x", "y", "z")])]), csv[[name]], tolerance = 1e-6)
    if (name %in% names(plotted_names)) {
      expect_equal(as.numeric(plotted[[plotted_names[[name]]]]), as.numeric(arr), tolerance = 1e-6)
    }
  }
  empty <- result
  empty$A <- empty$B <- 0L
  empty_summary <- ftq_exploratory_summary(empty, assignment_mode = "full_streamline", write_maps = FALSE)
  empty_html <- file.path(dirname(report$html), "no-targets.html")
  FibrTraQ:::ftq_write_exploratory_html_report(empty_html, empty_summary, list())
  expect_match(paste(readLines(empty_html), collapse = ""), "Not available (no eligible voxels)", fixed = TRUE)
})

test_that("exploratory categorical previews validate selected slices", {
  expect_error(FibrTraQ:::ftq_plot_categorical_montage(array(0, c(3, 3, 3)),
    data.frame(), tempfile(fileext = ".png"), "test", slices = 4), "slice|Slice")
  expect_error(FibrTraQ:::ftq_plot_categorical_montage(array(0, c(3, 3, 3)),
    data.frame(), tempfile(fileext = ".png"), "test", slices = 1.5), "slice|Slice")
})

test_that("target reports reuse current guides and keep complete target information", {
  root <- tempfile("target-guide-")
  figures <- file.path(root, "subject", "focal", "figures")
  dir.create(figures, recursive = TRUE)
  old_guide <- file.path(root, "OUTPUTS_README.html")
  writeLines("Existing customized guide", old_guide)
  guide <- FibrTraQ:::ftq_write_output_guide(root)
  expect_false(identical(guide, old_guide))
  input <- file.path(dirname(figures), "test_fibrtraq_per-voxel-hits.csv")
  metadata <- file.path(dirname(figures), "test_fibrtraq_run-metadata.json")
  utils::write.csv(public_exploratory_fixture(), input, row.names = FALSE)
  jsonlite::write_json(list(summary = list(target_intersection_mode = "full_streamline")), metadata)
  summary <- ftq_exploratory_summary(input, write_maps = FALSE)
  summary$target_totals <- summary$target_totals[rep(1, 26), ]
  summary$target_totals$target <- paste0("Target_", seq_len(26))
  output <- file.path(figures, "target-report.html")
  FibrTraQ:::ftq_write_exploratory_html_report(output, summary, list())
  html <- paste(readLines(output), collapse = "\n")
  expect_match(html, basename(guide), fixed = TRUE)
  expect_match(html, "../test_fibrtraq_run-metadata.json", fixed = TRUE)
  expect_match(html, "Target_26", fixed = TRUE)
  expect_false(grepl("download contains every target|Ranked preview", html))
  expect_identical(list.files(figures), "target-report.html")
  expect_identical(readLines(old_guide), "Existing customized guide")
})

test_that("standalone exploratory script forwards explicit method radius and slices", {
  skip_if_not_installed("RNifti")
  root <- tempfile("exploratory-script-")
  dir.create(root)
  input <- file.path(root, "result.csv")
  template <- file.path(root, "template.nii.gz")
  utils::write.csv(public_exploratory_fixture(), input, row.names = FALSE)
  RNifti::writeNifti(array(0, c(3, 3, 3)), template)
  script <- system.file("scripts", "fibrtraq-exploratory-report.R", package = "FibrTraQ")
  expect_true(nzchar(script))
  env <- new.env(parent = globalenv())
  args <- c(input, template, file.path(root, "report"), "--assignment-mode", "local_search",
    "--search-radius-mm", "3.25", "--slices", "1,3", "--prefix", "S1c")
  env$commandArgs <- function(trailingOnly = FALSE) args
  invisible(utils::capture.output(sys.source(script, envir = env)))
  expect_true(file.exists(env$report$html))
  expect_identical(env$report$summary$provenance$target_intersection_mode, "local_search")
  expect_identical(env$report$summary$provenance$local_search_radius_mm, 3.25)
  expect_equal(env$report_slices, c(1, 3))
  args <- c(input, template, file.path(root, "invalid"), "--search-radius-mm", "NaN")
  expect_error(sys.source(script, envir = env), "finite, positive number")
})
