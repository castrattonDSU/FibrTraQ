\name{ftq_stat_negative_control}
\alias{ftq_stat_negative_control}
\title{Compare correct runs against negative controls}
\usage{
ftq_stat_negative_control(
  correct,
  control,
  by,
  metric = "target_hit_rate",
  alternative = c("greater", "two.sided", "less")
)
}
\arguments{
\item{correct}{Data frame of subject/run summaries.}
\item{control}{Data frame of matched negative-control summaries.}
\item{by}{Column(s) used to pair rows.}
\item{metric}{Metric column to compare, or one of \code{"target_hit_rate"} and \code{"leave_no_target_fraction"} if raw count columns are present.}
\item{alternative}{Alternative hypothesis passed to \code{stats::wilcox.test()}.}
}
\value{
A list with paired rows, effect size, and paired Wilcoxon result.
}
\description{
Compare correct FibrTraQ runs against wrong-target, shifted-ROI, or other
negative-control runs in a paired design.
}
