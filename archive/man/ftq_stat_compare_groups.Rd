\name{ftq_stat_compare_groups}
\alias{ftq_stat_compare_groups}
\title{Compare subject-level metrics between groups}
\usage{
ftq_stat_compare_groups(
  data,
  metric,
  group_col = "group",
  subject_col = "subject",
  method = c("wilcox", "t")
)
}
\arguments{
\item{data}{Subject-level summary data.}
\item{metric}{Metric column to compare.}
\item{group_col}{Grouping column.}
\item{subject_col}{Subject column.}
\item{method}{\code{"wilcox"} or \code{"t"}.}
}
\value{
A data.frame with one row per metric comparison.
}
\description{
Run a simple two-group comparison on subject-level FibrTraQ metrics.
}
