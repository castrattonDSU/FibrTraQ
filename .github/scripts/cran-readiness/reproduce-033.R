# Optional diagnostic using the exact, immutable archive submitted as 0.3.3.
source(".github/scripts/cran-readiness/preflight.R")
evidence <- file.path(Sys.getenv("RUNNER_TEMP"), "ftq-cran033-path-reproduction")
dir.create(evidence, recursive = TRUE, showWarnings = FALSE)
archive <- file.path(Sys.getenv("RUNNER_TEMP"), "ftq-cran033-source", "FibrTraQ_0.3.3.tar.gz")
expected <- "579d0c0445b9e0b95705cc32de167d707381258576c4d1e1ec2fde0fb29a5713"
stopifnot(identical(digest::digest(archive, algo = "sha256", file = TRUE), expected))
writeLines(c(paste0("SHA256: ", expected),
             "SourceCommit: 8ba13990bdb1cc55730171427be44b4fb1fce082"),
           file.path(evidence, "baseline-source.txt"))
library <- file.path(evidence, "library")
dir.create(library, recursive = TRUE)
install_log <- file.path(evidence, "install.log")
status <- system2(file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", paste0("--library=", shQuote(library)), shQuote(archive)),
  stdout = install_log, stderr = install_log, timeout = 300)
if (!identical(status, 0L)) stop("The archived 0.3.3 package did not install; inspect install.log")
check_env <- ftq_cran_environment("long-temp", evidence)
result <- callr::r(function(library, evidence) {
  .libPaths(c(library, .libPaths()))
  library(FibrTraQ, lib.loc = library)
  stopifnot(identical(as.character(packageVersion("FibrTraQ", lib.loc = library)), "0.3.3"))
  output <- tempfile("fibrtraq-known-answer-test-")
  warnings <- character()
  error <- NULL
  answer <- withCallingHandlers(tryCatch(
    ftq_run_example(output_dir = output, backend = "cpp", verbose = FALSE),
    error = function(e) { error <<- conditionMessage(e); NULL }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  # Preserve partial files under short, flat names before R removes its tempdir.
  files <- list.files(output, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                      no.. = TRUE)
  files <- files[!dir.exists(files)]
  partial <- file.path(evidence, "partial-output")
  dir.create(partial)
  destinations <- file.path(partial, sprintf("%04d-%s", seq_along(files), basename(files)))
  copied <- if (length(files)) file.copy(files, destinations) else logical()
  manifest <- data.frame(Source = files, Characters = nchar(files),
                         PreservedAs = destinations, Copied = copied)
  utils::write.csv(manifest, file.path(evidence, "partial-output-files.csv"), row.names = FALSE)
  list(error = error, warnings = warnings, output_dir = output,
       tempdir = tempdir(), tempdir_characters = nchar(tempdir()),
       example_passed = isTRUE(answer$validation$passed), files = manifest)
}, args = list(library = library, evidence = evidence), env = check_env,
   user_profile = FALSE, system_profile = FALSE, timeout = 300,
   stdout = file.path(evidence, "example-stdout.txt"),
   stderr = file.path(evidence, "example-stderr.txt"))
jsonlite::write_json(result, file.path(evidence, "reproduction.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null")
writeLines(c(result$error, result$warnings), file.path(evidence, "captured-conditions.txt"))
diagnosis <- ftq_cran_report_failure(result)
jsonlite::write_json(diagnosis, file.path(evidence, "reproduction-disposition.json"),
                     auto_unbox = TRUE, pretty = TRUE)
if (!diagnosis$reproduced) {
  stop("Original failure was not reproduced with an identified >=260-character report path; inspect the diagnostic evidence.")
}
cat("Reproduced the original report-writing failure with an identified >=260-character path.\n")
