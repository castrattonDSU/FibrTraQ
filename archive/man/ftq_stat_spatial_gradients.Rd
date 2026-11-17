\name{ftq_stat_spatial_gradients}
\alias{ftq_stat_spatial_gradients}
\title{Estimate focal-voxel spatial gradients}
\usage{
ftq_stat_spatial_gradients(
  result_csv,
  target_labels = NULL,
  min_connected_voxels = 5L
)
}
\arguments{
\item{result_csv}{FibrTraQ per-voxel CSV path or data.frame.}
\item{target_labels}{Optional target columns to include.}
\item{min_connected_voxels}{Minimum nonzero voxels required per target.}
}
\value{
A data.frame of target-by-axis slopes.
}
\description{
Fit target-fraction gradients against focal voxel x/y/z coordinates. For
publication inference, run this per subject and test the resulting subject-level
slopes.
}
