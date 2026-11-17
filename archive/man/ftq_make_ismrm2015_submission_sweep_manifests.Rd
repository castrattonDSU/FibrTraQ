\name{ftq_make_ismrm2015_submission_sweep_manifests}
\alias{ftq_make_ismrm2015_submission_sweep_manifests}
\title{Generate ISMRM 2015 submission sweep manifests}
\usage{
ftq_make_ismrm2015_submission_sweep_manifests(
  scoring_dir,
  tractograms,
  output_dir,
  bundles = NULL,
  run_output_root = NULL,
  subject = "sub-ismrm2015",
  derivative_name = "fibrtraq-ismrm2015-submission-sweep",
  target_label = "tail",
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
\item{tractograms}{Character vector of submitted \code{.trk} files.}
\item{output_dir}{Directory where manifest YAML files should be written.}
\item{bundles}{Optional bundle names. Defaults to all configured bundles.}
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
\item{force}{Whether existing manifest files may be overwritten.}
}
\value{
A data.frame describing generated submission-sweep manifests.
}
\description{
Build manifests that run one or more submitted tractograms against every
configured ISMRM endpoint pair. This estimates which known bundle endpoints
are recovered by each submitted tractogram.
}
