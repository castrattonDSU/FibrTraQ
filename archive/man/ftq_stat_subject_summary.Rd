\name{ftq_stat_subject_summary}
\alias{ftq_stat_subject_summary}
\title{Summarize a FibrTraQ result for subject-level statistics}
\usage{
ftq_stat_subject_summary(
  result_csv,
  subject = NA_character_,
  group = NA_character_,
  roi = NA_character_,
  target_labels = NULL,
  endpoint_summary = NULL,
  endpoint_penetration = NULL
)
}
\arguments{
\item{result_csv}{FibrTraQ per-voxel CSV path or data.frame.}
\item{subject}{Optional subject identifier.}
\item{group}{Optional group/condition label.}
\item{roi}{Optional focal ROI label.}
\item{target_labels}{Optional target columns to include.}
\item{endpoint_summary}{Optional endpoint-summary CSV/data frame. When
supplied, an additive \code{endpoints} table is returned.}
\item{endpoint_penetration}{Optional endpoint-penetration CSV/data frame.
When supplied with \code{endpoint_summary}, an additive \code{penetration}
table is returned.}
}
\value{
A list with \code{roi}, \code{targets}, and \code{voxels} data.frames. When
endpoint inputs are supplied, the list also contains additive \code{endpoints}
and \code{penetration} members.
}
\description{
Collapse a per-voxel FibrTraQ table into ROI-level and target-level metrics
suitable for publication analyses. Use this before group comparisons so each
subject contributes one independent row per target.
}
