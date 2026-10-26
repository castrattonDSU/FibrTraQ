\name{ftq_null_shifted_rois}
\alias{ftq_null_shifted_rois}
\title{Write shifted ROI masks for spatial-null controls}
\usage{
ftq_null_shifted_rois(
  mask,
  output_dir,
  shifts,
  prefix = tools::file_path_sans_ext(basename(mask)),
  overwrite = FALSE
)
}
\arguments{
\item{mask}{Path to a NIfTI ROI mask.}
\item{output_dir}{Directory for shifted masks.}
\item{shifts}{Matrix/data.frame with columns \code{x}, \code{y}, \code{z}, or a list of numeric length-3 shifts in voxels.}
\item{prefix}{Output prefix.}
\item{overwrite}{Whether existing shifted masks may be replaced.}
}
\value{
A data.frame of written mask paths and shifts.
}
\description{
Create translated ROI masks for spatial-null and negative-control manifests.
Shifted voxels that leave the image bounds are dropped and new voxels are
zero-filled.
}
