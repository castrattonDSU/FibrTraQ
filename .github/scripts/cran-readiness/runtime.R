# CI-only completion and elapsed-time checks for the exact source archive.

ftq_cran_check_completion <- function(result, lines) {
  stage_ok <- function(stage) {
    pattern <- paste0("^\\* checking ", stage, " \\.\\.\\.( \\[[^]]+\\])? OK$")
    any(grepl(pattern, lines, perl = TRUE))
  }
  c(ProcessSucceeded = identical(result$status, 0L),
    NoTimeout = !isTRUE(result$timeout),
    CheckFinished = identical(tail(grep("^\\* ", lines, value = TRUE), 1L), "* DONE"),
    TestsPassed = stage_ok("tests"),
    PDFManualPassed = stage_ok("PDF version of manual"),
    HTMLManualPassed = stage_ok("HTML version of manual") &&
      !any(grepl("Skipping checking HTML validation", lines, fixed = TRUE)))
}

ftq_cran_check_runtime <- function(duration, elapsed, lines, evidence) {
  # This six-minute release budget includes every R CMD check stage. The longer
  # process timeout remains separate so a slow check can finish its diagnostics.
  budget <- 360
  valid_seconds <- function(value) {
    is.numeric(value) && !is.complex(value) && !is.object(value) &&
      length(value) == 1L && is.finite(value) && value >= 0
  }
  reported_valid <- valid_seconds(duration)
  measured_valid <- valid_seconds(elapsed)
  valid <- reported_valid && measured_valid
  checks <- c(ReportedDurationValid = reported_valid,
              MeasuredElapsedValid = measured_valid,
              # rcmdcheck starts its timer after archive staging. Allow one second
              # for the two clocks' resolution; neither can exceed the budget.
              DurationsConsistent = valid && duration <= elapsed + 1,
              WithinRuntimeBudget = valid && max(duration, elapsed) <= budget)
  record <- c(list(BudgetSeconds = budget,
                   ReportedDurationSeconds = if (reported_valid) duration else NA_real_,
                   MeasuredElapsedSeconds = if (measured_valid) elapsed else NA_real_),
              as.list(checks))
  write.dcf(as.data.frame(record), file.path(evidence, "check-runtime.dcf"))

  # Preserve detailed stage timings without conflating CPU time with elapsed time.
  pattern <- "^\\* checking (.+?) \\.\\.\\.(?: \\[([^]]+)\\])? (.*)$"
  matches <- regmatches(lines, regexec(pattern, lines, perl = TRUE))
  matches <- matches[lengths(matches) > 0L]
  stages <- data.frame(Stage = vapply(matches, `[`, "", 2L),
                       Timing = vapply(matches, `[`, "", 3L),
                       Result = vapply(matches, `[`, "", 4L))
  utils::write.csv(stages, file.path(evidence, "check-stage-timings.csv"), row.names = FALSE)
  invisible(checks)
}
