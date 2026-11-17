# Run the bundled 14-streamline example. In RStudio the report opens on completion.
local({
  inputs <- system.file("extdata", "fibrtraq-known-answer", "inputs",
                        package = "FibrTraQ", mustWork = TRUE)
  input <- function(name) file.path(inputs, paste0("sub-example_", name))

  demo <- FibrTraQ::fibrtraq_single(
    tractogram = input("tractogram.trk"),
    focal_mask = input("focal.nii.gz"),
    target_masks = c(cortex = input("target-cortex.nii.gz"),
                     brainstem = input("target-brainstem.nii.gz")),
    local_mask = input("local.nii.gz"),
    endpoint_contexts = c(analysis_domain = input("context.nii.gz")),
    reference_image = input("reference.nii.gz"),
    subject = "sub-example",
    focal_label = "focal",
    track_space = "voxel",
    axis_flip = c(FALSE, FALSE, FALSE),
    discretization = "round",
    output_dir = tempfile("fibrtraq-example-"),
    open = interactive()
  )
  stopifnot(isTRUE(demo$index$validation_passed))
  invisible(demo)
})
