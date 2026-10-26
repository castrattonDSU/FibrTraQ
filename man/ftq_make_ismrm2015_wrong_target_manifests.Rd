\name{ftq_make_ismrm2015_wrong_target_manifests}
\alias{ftq_make_ismrm2015_wrong_target_manifests}
\title{Generate ISMRM 2015 wrong-target control manifests}
\usage{
ftq_make_ismrm2015_wrong_target_manifests(
  scoring_dir,
  output_dir,
  bundles = NULL,
  target_bundles = NULL,
  max_targets_per_bundle = NULL,
  run_output_root = NULL,
  subject = "sub-ismrm2015",
  derivative_name = "fibrtraq-ismrm2015-wrong-targets",
  track_space = "auto",
  discretization = "any",
  write_maps = FALSE,
  write_report = FALSE,
  overwrite = FALSE,
  include_background = FALSE,
  combine_targets = FALSE,
  force = FALSE
)
}
\arguments{
\item{scoring_dir}{Path to the unpacked ISMRM 2015 scoring-data directory.}
\item{output_dir}{Directory where manifest YAML files should be written.}
\item{bundles}{Optional source bundle names. Defaults to all configured bundles.}
\item{target_bundles}{Optional bundle names whose tail endpoints should be used as wrong targets.}
\item{max_targets_per_bundle}{Optional integer cap for quick pilot runs.}
\item{run_output_root}{Optional root directory for FibrTraQ outputs.}
\item{subject}{Subject label written into each manifest.}
\item{derivative_name}{Derivative/output name written into each manifest.}
\item{track_space}{Track coordinate-space setting written into each manifest.}
\item{discretization}{Discretization setting written into each manifest.}
\item{write_maps}{Whether generated manifests should request NIfTI maps.}
\item{write_report}{Whether generated manifests should request HTML reports.}
\item{overwrite}{Whether generated manifests should permit overwriting core output files by default.}
\item{include_background}{Whether to include \code{t1.nii.gz} as \code{report.background_image}.}
\item{combine_targets}{If TRUE, generate one manifest per source bundle with all eligible wrong
targets included as separate target columns.}
\item{force}{Whether existing manifest files may be overwritten.}
}
\value{
A data.frame describing generated control manifests.
}
\description{
Build negative-control manifests by pairing each ground-truth bundle's
tractogram and head endpoint with other bundles' tail endpoints. These runs
test whether FibrTraQ preserves endpoint specificity rather than counting any
focal streamline as a target connection.
}
