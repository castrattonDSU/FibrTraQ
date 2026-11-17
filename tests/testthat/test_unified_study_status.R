test_that("the study report makes failed and unrecorded analyses visible at the top", {
  out <- tempfile("study-status-")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  path <- file.path(out, "report.html")
  runs <- data.frame(subject = c("first", "second"), focal_region = "focal",
                     run_status = c("completed", "failed"),
                     biological_qc_status = c("unavailable", "unavailable"))
  render <- function(rows, planned = nrow(rows)) {
    sections <- rep(list(list()), nrow(rows))
    ftq_write_unified_report(path, out, rows, sections, "test", planned_runs = planned)
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }

  partial <- render(runs, 4L)
  expect_match(partial, "Run records available for 2 of 4 planned analyses.", fixed = TRUE)
  expect_match(partial, "One analysis failed.", fixed = TRUE)
  expect_match(partial, "Select a run marked (failed)", fixed = TRUE)
  expect_match(partial, ">second / focal (failed)</option>", fixed = TRUE)
  expect_match(partial, ">first / focal</option>", fixed = TRUE)
  expect_true(regexpr("Study incomplete.", partial, fixed = TRUE)[1] <
              regexpr("<main id='ftq-active'>", partial, fixed = TRUE)[1])

  runs$run_status <- "completed"
  expect_false(grepl("Study incomplete.", render(runs), fixed = TRUE))
  expect_false(grepl("(failed)</option>", render(runs), fixed = TRUE))
  interrupted <- render(runs[1, , drop = FALSE], 2L)
  expect_match(interrupted, "Run records available for 1 of 2 planned analyses.", fixed = TRUE)
  expect_false(grepl("analysis failed.", interrupted, fixed = TRUE))
  empty <- render(runs[FALSE, , drop = FALSE], 2L)
  expect_match(empty, "Run records available for 0 of 2 planned analyses.", fixed = TRUE)
})
