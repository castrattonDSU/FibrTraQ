# Run from the repository root: Rscript .github/scripts/cran-readiness/test-runtime.R
# These CI-wrapper checks are excluded from the submitted source package.
source(".github/scripts/cran-readiness/runtime.R")

local({
  evidence <- tempfile("cran-runtime-guards-")
  dir.create(evidence)
  on.exit(unlink(evidence, recursive = TRUE))
  assertions <- 0L
  expect <- function(value, label) {
    assertions <<- assertions + 1L
    if (!isTRUE(value)) stop("Runtime guard regression: ", label)
  }
  lines <- c("* checking tests ... [210s/220s] OK",
             "* checking PDF version of manual ... [7s] OK",
             "* checking HTML version of manual ... [0s] OK", "* DONE")
  result <- list(status = 0L, timeout = FALSE)
  expect(all(ftq_cran_check_completion(result, lines)), "completed check")
  for (stage in c("tests", "PDF version of manual", "HTML version of manual")) {
    incomplete <- lines[!startsWith(lines, paste0("* checking ", stage, " ..."))]
    expect(!all(ftq_cran_check_completion(result, incomplete)), paste("missing", stage))
  }
  for (incomplete in list(character(), head(lines, -1L),
                          sub("tests ... [210s/220s] OK", "tests ... [210s/220s] ERROR", lines, fixed = TRUE),
                          c("Skipping checking HTML validation: no command tidy found", lines))) {
    expect(!all(ftq_cran_check_completion(result, incomplete)), "incomplete/failed/skipped check")
  }
  expect(!all(ftq_cran_check_completion(list(status = 1L), lines)), "failed process")
  expect(!all(ftq_cran_check_completion(list(status = 0L, timeout = TRUE), lines)), "timed-out process")
  expect(!all(ftq_cran_check_completion(list(), lines)), "missing process status")

  check <- function(duration, elapsed) ftq_cran_check_runtime(duration, elapsed, lines, evidence)
  for (value in c(0, 359.999, 360)) {
    expect(all(check(value, value)), paste("valid boundary", value))
  }
  expect(!all(check(360.001, 360.001)), "over six minutes")
  expect(!all(check(359, 360.001)), "measured time over budget")
  expect(!all(check(360.001, 360)), "reported time over budget")
  expect(!all(check(200, 100)), "inconsistent clocks")
  expect(all(check(100, 101)), "archive staging included in measured elapsed")
  malformed <- list(NULL, numeric(), c(1, 2), NA_real_, NaN, Inf, -Inf, -1,
                    "1", TRUE, list(1), 1 + 0i, as.difftime(1, units = "mins"))
  for (value in malformed) {
    expect(!all(check(value, 1)), "invalid reported duration")
    expect(!all(check(1, value)), "invalid measured elapsed")
  }
  check(359, 360)
  record <- read.dcf(file.path(evidence, "check-runtime.dcf"))
  expect(identical(unname(record[1L, c("BudgetSeconds", "ReportedDurationSeconds", "MeasuredElapsedSeconds")]),
                   c("360", "359", "360")), "persist exact seconds and budget")
  stages <- utils::read.csv(file.path(evidence, "check-stage-timings.csv"))
  expect(identical(stages$Timing, c("210s/220s", "7s", "0s")), "retain CPU/elapsed stage evidence")
  expect(identical(stages$Stage, c("tests", "PDF version of manual", "HTML version of manual")),
         "retain stage names")
  ftq_cran_check_runtime(1, 1, character(), evidence)
  expect(nrow(utils::read.csv(file.path(evidence, "check-stage-timings.csv"))) == 0L,
         "missing log remains empty evidence, not fabricated timing")

  # Exercise the real test-result structure and the shipped wrapper's CSV branch.
  fixture <- file.path(evidence, "fixture")
  dir.create(fixture)
  writeLines('testthat::test_that("timing fixture", testthat::expect_true(TRUE))',
             file.path(fixture, "test-timing.R"))
  results <- testthat::test_dir(fixture, reporter = "silent")
  wrapper <- normalizePath("tests/testthat.R")
  timing_file <- file.path(evidence, "test-timings.csv")
  previous <- Sys.getenv("FTQ_TEST_TIMINGS_FILE", unset = NA_character_)
  on.exit(if (is.na(previous)) Sys.unsetenv("FTQ_TEST_TIMINGS_FILE") else
    Sys.setenv(FTQ_TEST_TIMINGS_FILE = previous), add = TRUE)
  run_wrapper <- function(callback) {
    env <- new.env(parent = globalenv())
    env$library <- function(...) invisible(NULL)
    env$test_check <- callback
    sys.source(wrapper, envir = env)
  }
  Sys.unsetenv("FTQ_TEST_TIMINGS_FILE")
  runs <- 0L
  callback <- function(...) { runs <<- runs + 1L; results }
  run_wrapper(callback)
  expect(runs == 1L && !file.exists(timing_file), "ordinary run has no optional timing file")
  Sys.setenv(FTQ_TEST_TIMINGS_FILE = timing_file)
  run_wrapper(callback)
  timing <- utils::read.csv(timing_file)
  expect(runs == 2L && nrow(timing) == 1L && timing$passed == 1L && timing$real >= 0,
         "capture timing from exactly one execution")
  expect(identical(names(timing), c("file", "test", "user", "system", "real", "nb", "passed", "failed",
                                   "error", "warning", "skipped")), "stable timing columns")
  failed <- tryCatch(run_wrapper(function(...) stop("original test failure")), error = identity)
  expect(inherits(failed, "error") && conditionMessage(failed) == "original test failure",
         "preserve test failure behavior")
  cat("CRAN runtime guard assertions passed:", assertions, "\n")
})
