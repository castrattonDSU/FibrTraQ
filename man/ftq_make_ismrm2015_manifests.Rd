\name{ftq_make_ismrm2015_manifests}
\alias{ftq_make_ismrm2015_manifests}
\title{Generate ISMRM 2015 validation manifests}
\usage{
ftq_make_ismrm2015_manifests(
  scoring_dir,
  output_dir,
  bundles = NULL,
  run_output_root = NULL,
  subject = "sub-ismrm2015",
  derivative_name = "fibrtraq-ismrm2015",
  target_label = "tail",
  track_space = "auto",
  discretization = "any",
  write_maps = TRUE,
  write_report = TRUE,
  overwrite = FALSE,
  include_background = TRUE,
  force = FALSE
)
}
\arguments{
\item{scoring_dir}{Path to the unpacked ISMRM 2015 scoring-data directory.
This directory should contain \code{config_file_segmentation.json}, \code{bundles/},
\code{ROI/}, and usually \code{t1.nii.gz}.}

\item{output_dir}{Directory where manifest YAML files should be written.}

\item{bundles}{Optional character vector of bundle names to generate. If
omitted, all configured bundles with matching \code{.trk} files are generated.}

\item{run_output_root}{Optional root directory for FibrTraQ outputs referenced
inside the generated manifests. Defaults to
\code{file.path(scoring_dir, "derivatives", derivative_name)}.}

\item{subject}{Subject label written into each manifest.}

\item{derivative_name}{Derivative/output name written into each manifest.}

\item{target_label}{Target ROI label used in the per-bundle manifest.}

\item{track_space}{Track coordinate-space setting written into each manifest.}

\item{discretization}{Discretization setting written into each manifest.}

\item{write_maps}{Whether generated manifests should request NIfTI maps.}

\item{write_report}{Whether generated manifests should request HTML reports.}

\item{overwrite}{Whether generated manifests should permit overwriting core
output files by default.}

\item{include_background}{Whether to include \code{t1.nii.gz} as
\code{report.background_image} when it exists.}

\item{force}{Whether existing manifest files may be overwritten.}
}
\value{
A data.frame describing generated manifests and input paths.
}
\description{
Read the ISMRM 2015 scoring-data segmentation configuration and create one
FibrTraQ manifest per bundle/sub-bundle with an available ground-truth \code{.trk}
file and configured head/tail endpoint masks.
}
