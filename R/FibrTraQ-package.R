#' FibrTraQ
#'
#' Voxel-resolved focal-region target profiling from TrackVis (`.trk`) or
#' MRtrix (`.tck`) diffusion tractograms and aligned NIfTI masks. Both formats
#' support the same summaries, density maps, HTML reports and interactive
#' three-dimensional quality-control viewer. Spatial alignment must be
#' independently verified before interpreting results.
#'
#' @seealso [ftq_run_example()] for a bundled known-answer example and
#'   [ftq_run_files()] for analysis of existing input files.
#'
#' @keywords internal
#' @useDynLib FibrTraQ, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
