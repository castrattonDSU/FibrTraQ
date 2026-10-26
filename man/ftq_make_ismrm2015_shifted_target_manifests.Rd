\name{ftq_make_ismrm2015_shifted_target_manifests}
\alias{ftq_make_ismrm2015_shifted_target_manifests}
\title{Generate ISMRM 2015 shifted-target control manifests}
\usage{
ftq_make_ismrm2015_shifted_target_manifests(
  scoring_dir,
  output_dir,
  shifted_roi_dir,
  shift = c(8L, 0L, 0L),
  bundles = NULL,
  run_output_root = NULL,
  subject = "sub-ismrm2015",
  derivative_name = "fibrtraq-ismrm2015-shifted-targets",
  target_label = "shifted_tail",
  track_space = "auto",
  discretization = "any",
  write_maps = FALSE,
  write_report = FALSE,
  overwrite = FALSE,
  include_background = FALSE,
  force = FALSE
)
}
\arguments{
\item{scoring_dir}{Path to the unpacked ISMRM 2015 scoring-data directory.}
\item{output_dir}{Directory where manifest YAML files should be written.}
\item{shifted_roi_dir}{Directory where shifted target masks should be written.}
\item{shift}{Integer length-3 vector \code{c(dx, dy, dz)} in voxel units.}
\item{bundles}{Optional source bundle names. Defaults to all configured bundles.}
\item{run_output_root}{Optional root directory for FibrTraQ outputs.}
\item{subject}{Subject label written into each manifest.}
\item{derivative_name}{Derivative/output name written into each manifest.}
\item{target_label}{Target ROI label used in generated manifests.}
\item{track_space}{Track coordinate-space setting written into each manifest.}
\item{discretization}{Discretization setting written into each manifest.}
\item{write_maps}{Whether generated manifests should request NIfTI maps.}
\item{write_report}{Whether generated manifests should request HTML reports.}
\item{overwrite}{Whether generated manifests should permit overwriting core output files by default.}
\item{include_background}{Whether to include \code{t1.nii.gz} as \code{report.background_image}.}
\item{force}{Whether existing manifest and shifted ROI files may be overwritten.}
}
\value{
A data.frame describing generated shifted-target manifests.
}
\description{
Create non-wrapping shifted copies of each configured target endpoint mask and
build one control manifest per bundle. These controls test spatial
specificity: moving the target ROI should strongly reduce target associations.
}
