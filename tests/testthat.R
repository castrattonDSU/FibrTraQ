library(testthat)
library(FibrTraQ)

results <- test_check("FibrTraQ")

# Optional CI evidence from this same run; the normal check reporter is unchanged.
timing_file <- Sys.getenv("FTQ_TEST_TIMINGS_FILE", "")
if (nzchar(timing_file)) {
  timings <- as.data.frame(results)
  columns <- c("file", "test", "user", "system", "real", "nb", "passed", "failed",
               "error", "warning", "skipped")
  utils::write.csv(timings[, columns, drop = FALSE], timing_file, row.names = FALSE)
}
