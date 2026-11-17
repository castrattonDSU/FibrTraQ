# R CMD check runs this in a fresh installed-package session.
library(FibrTraQ)
commands <- c("fibrtraq_single", "fibrtraq_study")
stopifnot(identical(sort(getNamespaceExports("FibrTraQ")), commands),
          identical(sort(ls("package:FibrTraQ")), commands))

inputs <- system.file("extdata", "fibrtraq-known-answer", "inputs",
                      package = "FibrTraQ", mustWork = TRUE)
input <- function(name) file.path(inputs, paste0("sub-example_", name))
tractogram <- input("tractogram.trk")
masks <- c(focal = input("focal.nii.gz"), cortex = input("target-cortex.nii.gz"),
           brainstem = input("target-brainstem.nii.gz"))
single_root <- tempfile("s-")
study_root <- tempfile("g-")
single <- fibrtraq_single(
  tractogram = tractogram, focal_mask = masks[["focal"]], target_masks = masks[-1L],
  subject = "s1", focal_label = "focal", local_mask = "union",
  output_dir = single_root, reference_image = input("reference.nii.gz"),
  track_space = "voxel", discretization = "round", verbose = FALSE, progress = FALSE
)
# Two labels share the synthetic fixture so both public calls have the same expected counts.
study <- fibrtraq_study(
  tractograms = c(s1 = tractogram, s2 = tractogram), region_masks = list(s1 = masks, s2 = masks),
  focal_regions = "focal", local_masks = "union",
  reference_images = c(s1 = input("reference.nii.gz"), s2 = input("reference.nii.gz")),
  output_dir = study_root, track_space = "voxel", discretization = "round",
  verbose = FALSE, progress = FALSE
)

expected_union <- Reduce(`|`, lapply(masks, function(path) RNifti::readNifti(path) != 0))
for (run in list(single, study)) {
  root <- run$output_dir
  pages <- list.files(root, "[.]html$", recursive = TRUE, full.names = TRUE)
  scenes <- list.files(root, "_qc-scene[.]js$", recursive = TRUE)
  stopifnot(all(run$index$validation_passed), all(run$index$run_status == "completed"),
            length(pages) == 1L, identical(normalizePath(pages), normalizePath(run$report_index)),
            length(scenes) == nrow(run$index),
            length(run$prepared_inputs$local_mask) == length(unique(run$index$subject)),
            identical(normalizePath(dirname(run$index$result_csv)),
                      normalizePath(file.path(root, run$index$subject, run$index$focal_region, "tables"))))
  for (path in run$prepared_inputs$local_mask) {
    actual_union <- RNifti::readNifti(path)
    stopifnot(identical(dim(actual_union), dim(expected_union)), all(actual_union == expected_union))
  }
  html <- paste(readLines(pages, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  stopifnot(grepl("Open 3D viewer", html, fixed = TRUE))
  for (scene in scenes) stopifnot(file.exists(file.path(root, scene)))
}
stopifnot(nrow(study$index) == 2L, nrow(single$index) == 1L)
for (field in c("result_csv", "rsc_csv", "fate_csv")) {
  first <- utils::read.csv(single$index[[field]])
  for (path in study$index[[field]]) stopifnot(identical(first, utils::read.csv(path)))
}
