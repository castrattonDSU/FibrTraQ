library(testthat)
library(FibrTraQ)

ftq_report_link_test_index <- function(reports) {
  data.frame(
    subject = "subject1", focal_region = "region1", run_status = "completed",
    validation_passed = TRUE, biological_qc_status = "not_reviewed",
    primary_report = reports[[1L]], exploratory_report = reports[[2L]],
    interactive_3d_qc = reports[[3L]], stringsAsFactors = FALSE
  )
}

ftq_report_link_test_hrefs <- function(html) {
  links <- regmatches(html, gregexpr("href=['\"][^'\"]+['\"]", html))[[1L]]
  sub("['\"]$", "", sub("^href=['\"]", "", links))
}

test_that("URL path components preserve separators and encode raw UTF-8 once", {
  expect_identical(
    FibrTraQ:::ftq_url_path("one two/#%/%20/caf\u00e9?.html"),
    "one%20two/%23%25/%2520/caf%C3%A9%3F.html"
  )
  expect_identical(
    FibrTraQ:::ftq_url_path("a&b/'quoted'\"+/-._~"),
    "a%26b/%27quoted%27%22%2B/-._~"
  )
  expect_identical(FibrTraQ:::ftq_url_path("../plots/a b.svg"),
                   "../plots/a%20b.svg")
  expect_identical(FibrTraQ:::ftq_url_path("/a/%23/"), "/a/%2523/")
})

test_that("Windows drive paths are handled lexically on every host", {
  roots <- c("C:\\study\\output", "C:/study/output", "c:/study/output/")
  paths <- c("C:\\study\\output\\subject 1\\report#%20.html",
             "C:/study/output/subject 1/report#%20.html")
  for (root in roots) for (path in paths) {
    expect_identical(FibrTraQ:::ftq_report_relative_path(path, root),
                     "subject 1/report#%20.html")
    expect_identical(FibrTraQ:::ftq_report_href(path, root),
                     "subject%201/report%23%2520.html")
  }
  expect_identical(
    FibrTraQ:::ftq_report_relative_path("C:\\study\\reports\\a.html",
                                        "C:\\study\\indexes\\nested\\"),
    "../../reports/a.html"
  )
  expect_identical(FibrTraQ:::ftq_report_relative_path("C:/study/output", "C:/study/output/"),
                   ".")
  expect_identical(FibrTraQ:::ftq_report_relative_path("D:\\reports\\a b.html", "C:/study"),
                   "D:/reports/a b.html")
  expect_identical(FibrTraQ:::ftq_report_href("D:\\reports\\a b.html", "C:/study"),
                   "file:///D:/reports/a%20b.html")
})

test_that("UNC paths stay relative only within the same server and share", {
  for (root in c("\\\\server\\share\\study\\", "//server/share/study/")) {
    expect_identical(
      FibrTraQ:::ftq_report_relative_path("\\\\server\\share\\study\\reports\\a.html", root),
      "reports/a.html"
    )
    expect_identical(FibrTraQ:::ftq_report_href("//server/share/other/a #.html", root),
                     "../other/a%20%23.html")
    expect_identical(FibrTraQ:::ftq_report_relative_path("//server/other/a.html", root),
                     "//server/other/a.html")
    expect_identical(FibrTraQ:::ftq_report_href("//server/other/a #.html", root),
                     "file://server/other/a%20%23.html")
    expect_identical(FibrTraQ:::ftq_report_href("//elsewhere/share/a.html", root),
                     "file://elsewhere/share/a.html")
  }
})

test_that("native Windows normalization cannot reintroduce UNC backslashes", {
  normalized_inputs <- list()
  windows_platform <- .Platform
  windows_platform$OS.type <- "windows"
  normalize <- FibrTraQ:::ftq_report_normalize_path
  environment(normalize) <- list2env(list(
    .Platform = windows_platform,
    normalizePath = function(path, winslash, mustWork) {
      normalized_inputs[[length(normalized_inputs) + 1L]] <<-
        list(path = path, winslash = winslash, mustWork = mustWork)
      # Observed on Windows even when winslash = "/" was requested.
      paste0("\\\\", substring(path, 3L))
    }
  ), parent = environment(normalize))
  for (path in c("\\\\server\\share\\report %20.html",
                 "//server/share/report %20.html")) {
    expect_identical(normalize(path), "//server/share/report %20.html")
  }
  expect_identical(normalized_inputs, rep(list(list(
    path = "//server/share/report %20.html", winslash = "/", mustWork = FALSE
  )), 2L))
})

test_that("file URIs distinguish drive, UNC, and native absolute paths", {
  expect_identical(FibrTraQ:::ftq_file_url("C:\\study space\\a#%20.html"),
                   "file:///C:/study%20space/a%23%2520.html")
  expect_identical(FibrTraQ:::ftq_file_url("C:/study space/a#%20.html"),
                   "file:///C:/study%20space/a%23%2520.html")
  expect_identical(FibrTraQ:::ftq_file_url("\\\\server\\share\\caf\u00e9.html"),
                   "file://server/share/caf%C3%A9.html")
  if (.Platform$OS.type != "windows") {
    expect_identical(FibrTraQ:::ftq_file_url("/study space/a#%20.html"),
                     "file:///study%20space/a%23%2520.html")
    expect_identical(FibrTraQ:::ftq_report_relative_path("/runner/a.html", "/run/"),
                     "../runner/a.html")
    expect_identical(FibrTraQ:::ftq_report_relative_path("/reports/a.html", "/"),
                     "reports/a.html")
  }
})

test_that("a literal POSIX backslash is encoded as a filename character", {
  expect_identical(FibrTraQ:::ftq_url_path("report\\name.html"), "report%5Cname.html")
  if (.Platform$OS.type != "windows") {
    root <- tempfile("report-backslash-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    path <- file.path(root, "report\\name.html")
    writeLines("literal backslash", path)
    expect_identical(FibrTraQ:::ftq_report_relative_path(path, root), "report\\name.html")
    href <- FibrTraQ:::ftq_report_href(path, root)
    expect_identical(href, "report%5Cname.html")
    expect_identical(readLines(file.path(root, utils::URLdecode(href))), "literal backslash")
  }
})

test_that("regional index handles native Windows strings and the original schema", {
  # Capture only the write boundary; all path processing uses the real writer.
  html <- NULL
  writer <- FibrTraQ:::ftq_write_regional_study_index
  environment(writer) <- list2env(list(
    writeLines = function(text, con, ...) {
      html <<- paste(text, collapse = "\n")
    },
    # The synthetic Windows paths are lexical fixtures, not writable locations.
    ftq_write_output_guide = function(output_dir) {
      paste0(output_dir, "/OUTPUTS_README.html")
    }
  ), parent = environment(writer))
  reports <- c("C:\\study #\\output\\subject 1\\standard %20.html",
               "C:\\study #\\output\\subject 1\\exploratory.html",
               "C:\\study #\\output\\subject 1\\QC.html")
  for (native_windows in c(TRUE, FALSE)) {
    index_path <- "C:\\study #\\output\\indexes\\index.html"
    output_dir <- "C:\\study #\\output"
    report_paths <- reports
    if (!native_windows) {
      index_path <- gsub("\\", "/", index_path, fixed = TRUE)
      output_dir <- gsub("\\", "/", output_dir, fixed = TRUE)
      report_paths <- gsub("\\", "/", report_paths, fixed = TRUE)
    }
    writer(ftq_report_link_test_index(report_paths), path = index_path,
           output_dir = output_dir, description = "link-regression")
    expect_identical(ftq_report_link_test_hrefs(html),
                     c("../subject%201/standard%20%2520.html",
                       "../subject%201/exploratory.html", "../subject%201/QC.html",
                       "../OUTPUTS_README.html"))
    expect_match(html, "not recorded", fixed = TRUE)
  }
})

test_that("raw regional hrefs survive relocation of the complete report tree", {
  parent <- tempfile("report-relocation-")
  original <- file.path(parent, "original tree # %20")
  relocated <- file.path(parent, "relocated tree # %20")
  dir.create(file.path(original, "indexes"), recursive = TRUE)
  dir.create(file.path(original, "subject 1 #"))
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  names <- c("standard %20.html", "exploratory #.html", "QC caf\u00e9.html")
  reports <- file.path(original, "subject 1 #", names)
  for (i in seq_along(reports)) writeLines(paste("report", i), reports[[i]])
  index <- file.path(original, "indexes", "index.html")
  FibrTraQ:::ftq_write_regional_study_index(
    ftq_report_link_test_index(reports), path = index, output_dir = original,
    description = "relocation-regression"
  )
  html <- paste(readLines(index, warn = FALSE), collapse = "\n")
  hrefs <- ftq_report_link_test_hrefs(html)
  expect_identical(hrefs, c("../subject%201%20%23/standard%20%2520.html",
                           "../subject%201%20%23/exploratory%20%23.html",
                           "../subject%201%20%23/QC%20caf%C3%A9.html",
                           "../OUTPUTS_README.html"))
  guide_content <- readLines(file.path(original, "OUTPUTS_README.html"), warn = FALSE)
  expect_match(paste(guide_content, collapse = "\n"), "FibrTraQ output guide", fixed = TRUE)
  expect_true(file.rename(original, relocated))
  expect_false(dir.exists(original))
  moved_index <- file.path(relocated, "indexes", "index.html")
  expect_identical(paste(readLines(moved_index, warn = FALSE), collapse = "\n"), html)
  for (i in seq_along(hrefs)) {
    target <- file.path(dirname(moved_index), utils::URLdecode(hrefs[[i]]))
    expect_true(file.exists(target))
    expect_identical(readLines(target, warn = FALSE),
                     if (i <= length(reports)) paste("report", i) else guide_content)
    expect_true(startsWith(normalizePath(target, winslash = "/"),
                           paste0(normalizePath(relocated, winslash = "/"), "/")))
  }
})
