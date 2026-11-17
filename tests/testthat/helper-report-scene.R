ftq_test_viewer_scene <- function(focal_dir) {
  files <- list.files(file.path(focal_dir, "figures"), "_qc-scene[.]js$", full.names = TRUE)
  expect_length(files, 1L)
  content <- readLines(files, warn = FALSE, encoding = "UTF-8")
  expect_length(content, 1L)
  prefix <- "document.currentScript.ftqScene = "
  expect_true(startsWith(content, prefix) && endsWith(content, ";"))
  jsonlite::fromJSON(substr(content, nchar(prefix) + 1L, nchar(content) - 1L),
                     simplifyVector = FALSE)
}
