#' FibrTraQ
#'
#' Voxel-resolved focal-region target profiling from TrackVis (`.trk`) or
#' MRtrix (`.tck`) diffusion tractograms and aligned NIfTI masks. Both formats
#' support the same summaries, density maps, HTML reports and interactive
#' three-dimensional quality-control viewer. Spatial alignment must be
#' independently verified before interpreting results.
#'
#' @seealso [fibrtraq_single()] for one subject and [fibrtraq_study()] for a study.
#'
#' @keywords internal
#' @useDynLib FibrTraQ, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
