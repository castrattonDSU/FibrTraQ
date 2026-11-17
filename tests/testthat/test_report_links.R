library(testthat)
library(FibrTraQ)

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
