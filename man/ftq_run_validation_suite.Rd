\name{ftq_run_validation_suite}
\alias{ftq_run_validation_suite}
\title{Run a suite of FibrTraQ validation manifests}
\usage{
ftq_run_validation_suite(
  manifests,
  output_csv = NULL,
  preflight_only = FALSE,
  overwrite = FALSE,
  write_maps = NULL,
  write_report = NULL,
  continue_on_error = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{manifests}{Character vector of manifest paths, or a directory containing
\code{.yml}, \code{.yaml}, or \code{.json} manifests. Directories are searched
recursively.}

\item{output_csv}{Optional path where the aggregate suite table should be
written.}

\item{preflight_only}{If TRUE, only run manifest preflight checks.}

\item{overwrite}{Whether existing FibrTraQ outputs may be overwritten during
full runs.}

\item{write_maps}{Optional override for manifest map output settings.}

\item{write_report}{Optional override for manifest report output settings.}

\item{continue_on_error}{If TRUE, failed manifests are retained as rows in the
output table. If FALSE, the first failed manifest stops the suite.}

\item{verbose}{Emit progress messages from the underlying FibrTraQ run.}
}
\value{
A data.frame with one row per manifest and summary columns for preflight
status, validation status, core output paths, and run totals.
}
\description{
Run preflight checks, and optionally full FibrTraQ analyses, across one or
more manifest files. This is intended for publication validation where each
dataset/application is represented by a manifest and the suite output becomes
an auditable summary table.
}
