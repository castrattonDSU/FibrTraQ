# Run from the repository root with the matching FibrTraQ package installed.
library(testthat)
archive_root <- normalizePath("archive", mustWork = TRUE)
Sys.setenv(FIBRTRAQ_REPO_ROOT = dirname(archive_root))
research <- source(file.path(archive_root, "load.R"))$value
test_env <- new.env(parent = research)
test_env$archive_root <- archive_root
test_env$research <- research
for (helper in list.files("tests/testthat", "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(helper, envir = test_env)
}
test_dir(file.path(archive_root, "tests", "testthat"), env = test_env,
         reporter = "summary", stop_on_failure = TRUE, stop_on_warning = TRUE)
