test_that("a study sheet preserves identifiers and resolves paths from the sheet", {
  root <- tempfile("sheet-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  root <- normalizePath(root)
  sheet <- file.path(root, "subjects.csv")
  absolute <- file.path(root, "target.nii.gz")
  file.create(absolute)
  writeLines(c(
    "\ufeffsubject,tractogram,reference_image,mask:focal,mask:target",
    paste0('001,"input/a,1.trk",,input\\focal.nii.gz,"', absolute, '"'),
    "NA,input/b.tck,image.nii.gz,input/focal.nii.gz,target.nii.gz"
  ), sheet, useBytes = TRUE)
  x <- ftq_read_study_sheet(sheet)
  expect_identical(names(x$tractograms), c("001", "NA"))
  expect_identical(basename(x$tractograms[["001"]]), "a,1.trk")
  normalized <- function(path) chartr("\\", "/", path)
  expect_equal(normalized(dirname(x$tractograms[["NA"]])), normalized(file.path(root, "input")))
  expect_identical(normalized(x$region_masks[["001"]][["focal"]]),
                   normalized(file.path(root, "input", "focal.nii.gz")))
  expect_identical(x$region_masks[["001"]][["target"]], normalizePath(absolute))
  expect_identical(names(x$reference_images), c("001", "NA"))
  expect_null(x$reference_images[["001"]])
  expect_identical(normalized(x$reference_images[["NA"]]), normalized(file.path(root, "image.nii.gz")))

  example <- system.file("examples", "study-sheet.csv", package = "FibrTraQ", mustWork = TRUE)
  expect_identical(names(ftq_read_study_sheet(example)$region_masks), c("001", "002"))
})

test_that("invalid study sheet columns and cells identify the problem", {
  sheet <- tempfile(fileext = ".csv")
  on.exit(unlink(sheet), add = TRUE)
  read <- function(header, rows = "001,one.trk,a.nii,b.nii") {
    writeLines(c(header, rows), sheet)
    ftq_read_study_sheet(sheet)
  }
  header <- "subject,tractogram,mask:A,mask:B"
  expect_error(ftq_read_study_sheet(sheet), "existing CSV")
  expect_error(read("id,tractogram,mask:A,mask:B"), "missing required columns: subject", fixed = TRUE)
  expect_error(read("subject,trk,mask:A,mask:B"), "missing required columns: tractogram", fixed = TRUE)
  expect_error(read("subject,tractogram,mask:A,unexpected"), "Unrecognized.*unexpected")
  expect_error(read("subject,tractogram,mask:A,mask:A"), "unique header")
  expect_error(read("subject,tractogram,mask:,mask:B"), "two named region columns")
  expect_error(read("subject,tractogram,mask:A,mask:a"), "distinct")
  expect_error(read(header, character()), "no subject rows")
  expect_error(read(header, c("001,t,a,b", "001,u,a,b")), "Repeated subject IDs: 001")
  expect_error(read(header, c("S1,t,a,b", "s1,u,a,b")), "distinct")
  expect_error(read(header, "001,t,,b"), "row 2 (mask:A)", fixed = TRUE)
  expect_error(read(header, c("001,t,a,b", "002,,a,b")), "row 3 (tractogram)", fixed = TRUE)
  expect_error(read(header, "001,t,a"), "Could not read the study sheet")
  expect_null(read(header)$reference_images)
})

test_that("check-only identifies missing sheet files for every affected subject", {
  root <- tempfile("sh-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  sheet <- file.path(root, "inputs.csv")
  writeLines(c("subject,tractogram,mask:A,mask:B",
               "001,missing.trk,a.nii,b.nii", "002,missing.tck,c.nii,d.nii"), sheet)
  output <- file.path(root, "results")
  x <- fibrtraq_study(study_sheet = sheet, focal_regions = "A", output_dir = output,
                      check_only = TRUE, verbose = FALSE)
  expect_false(x$passed)
  expect_identical(x$index$subject, c("001", "002"))
  expect_false(any(x$index$passed))
  expect_true(any(grepl("missing.trk", x$checks$detail[x$checks$subject == "001"], fixed = TRUE)))
  expect_true(any(grepl("missing.tck", x$checks$detail[x$checks$subject == "002"], fixed = TRUE)))
  expect_false(file.exists(output))
})

test_that("mixed-format CSV and named-list studies return the same results", {
  root <- tempfile("sh-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir <- file.path(root, "data")
  files <- tck_fixture_files(dir)
  sheet <- file.path(root, "study.csv")
  writeLines(c(
    "subject,tractogram,reference_image,mask:focal,mask:target_a,mask:target_b",
    "001,data/tracks.trk,data/reference.nii.gz,data/focal.nii.gz,data/target_a.nii.gz,data/target_b.nii.gz",
    "002,data/tracks.tck,,data/focal.nii.gz,data/target_a.nii.gz,data/target_b.nii.gz"
  ), sheet)
  settings <- list(focal_regions = "focal", discretization = "round", verbose = FALSE,
                   endpoint_analysis = "none", write_report = FALSE,
                   write_exploratory_report = FALSE, write_qc_viewer = FALSE)
  run <- function(input, destination) {
    do.call(fibrtraq_study, c(input, settings, list(output_dir = file.path(root, destination))))
  }
  from_sheet <- run(list(study_sheet = sheet), "s")
  masks <- c(focal = files$focal, files$targets)
  lists <- list(tractograms = c("001" = files$trk, "002" = files$tck),
                region_masks = list("001" = masks, "002" = masks),
                reference_images = list("001" = files$reference, "002" = NULL))
  from_lists <- run(lists, "l")
  expect_identical(from_sheet$index$subject, c("001", "002"))
  expect_true(all(from_sheet$index$validation_passed))
  expect_identical(from_sheet$index$run_status, from_lists$index$run_status)
  expect_true(file.exists(from_sheet$report_index))
  for (field in c("result_csv", "rsc_csv", "fate_csv")) {
    for (i in seq_len(2L)) {
      expect_identical(readLines(from_sheet$index[[field]][i]), readLines(from_lists$index[[field]][i]))
    }
  }
  for (conflict in names(lists)) {
    expect_error(run(c(list(study_sheet = sheet), lists[conflict]), "conflict"), "study_sheet")
  }
  expect_false(dir.exists(file.path(root, "conflict")))
})
