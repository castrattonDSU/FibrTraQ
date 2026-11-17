#' FibrTraQ (branded entry point)
#'
#' @description Branded wrapper around \code{\link{FibrTraQ}}.
#' @param ... Additional arguments forwarded to \code{\link{FibrTraQ}}.
#'   This is primarily for backward compatibility. Unknown arguments will error
#'   in \code{FibrTraQ()} because it does not accept arbitrary \code{...}.
#' @return Per-voxel table of focal-voxel target associations.
#' @noRd
fibrtraq <- function(...) FibrTraQ(...)

#' Summarize voxel-resolved target associations and RSC
#'
#' @description Wrapper around \code{\link{summary_connectivity}}.
#' @param ... Additional arguments forwarded to \code{\link{summary_connectivity}}.
#'   Unknown arguments will error in \code{summary_connectivity()} because it does not
#'   accept arbitrary \code{...}.
#' @return See \code{\link{summary_connectivity}}.
#' @noRd
rsc_summary <- function(...) summary_connectivity(...)
