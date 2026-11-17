#' Build tractogram-derived spatial context for reports
#'
#' Stream a TrackVis `.trk` or MRtrix `.tck` file and build a streamline-occupancy
#' density image on a reference voxel grid. This is a visualization and
#' quality-control context only: it is not an anatomical image and must not be
#' interpreted as brain, tissue, or lesion volume.
#'
#' Each streamline contributes at most one count to a voxel under the selected
#' stored-point or segment-interpolation sampling and discretization contract.
#'
#' @param tractogram TrackVis `.trk` or MRtrix `.tck` file.
#' @param template_image NIfTI/MGH image defining the output voxel grid.
#' @param output_dir Directory where a `context/` subdirectory is written. If
#'   omitted, the template image directory is used.
#' @param prefix Output filename prefix.
#' @param track_space Track coordinate space. Use the same setting as the
#'   FibrTraQ run whenever possible. `"voxelmm"` uses TrackVis header voxel
#'   sizes and the TrackVis corner-to-center convention. Native `.tck` files
#'   use RAS+ world millimeters; only `"auto"` or `"world"` with no axis flips
#'   is accepted, and the template/reference must declare valid spatial geometry.
#' @param axis_flip Logical length-3 vector applied after discretization.
#' @param track_chunk_size Number of streamlines read per chunk.
#' @param density_discretization Discretization used for the visualization
#'   density. This does not affect FibrTraQ counting outputs.
#' @param streamline_voxelization_mode Whether density uses stored sample
#'   points only or interpolates segments at a maximum 0.5-voxel axis step.
#' @param ref_image Optional reference image for world/TrackVis conversion.
#'   Defaults to `template_image`.
#' @param overwrite Whether existing context outputs may be replaced.
#'
#' @return A list with `density_image`, `summary_csv`, `metadata_json`,
#'   `summary`, and `metadata`.
#' @export
ftq_tractogram_context <- function(tractogram,
                                   template_image,
                                   output_dir = NULL,
                                   prefix = "fibrtraq",
                                   track_space = c("auto", "voxelmm", "vox10", "voxel", "world", "trackvis"),
                                   axis_flip = c(FALSE, FALSE, FALSE),
                                   track_chunk_size = 50000L,
                                   density_discretization = c("round", "floor", "ceil"),
                                   streamline_voxelization_mode = c("sample_points", "segment_interpolation"),
                                   ref_image = NULL,
                                   overwrite = FALSE) {
  ftq_build_tractogram_context(
    tractogram = tractogram,
    template_image = template_image,
    output_dir = output_dir,
    prefix = prefix,
    track_space = track_space,
    axis_flip = axis_flip,
    track_chunk_size = track_chunk_size,
    density_discretization = density_discretization,
    streamline_voxelization_mode = streamline_voxelization_mode,
    ref_image = ref_image,
    overwrite = overwrite
  )
}
