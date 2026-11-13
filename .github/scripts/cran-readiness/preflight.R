# CI-only checks of the exact source archive and CRAN checking environment.

ftq_cran_readme_links <- function(source, evidence) {
  readme <- file.path(source, "README.md")
  stopifnot(file.exists(readme))
  pandoc <- unname(Sys.which("pandoc"))
  if (!nzchar(pandoc)) stop("Pandoc is required for the archive README link check")
  parsed <- tempfile(fileext = ".json")
  on.exit(unlink(parsed))
  status <- system2(pandoc, c("--from=gfm", "--to=json", shQuote(readme)),
                    stdout = parsed)
  if (!identical(status, 0L)) stop("Could not parse the archived README")
  document <- jsonlite::fromJSON(parsed, simplifyVector = FALSE)
  links <- character()
  visit <- function(node) {
    if (!is.list(node)) return(invisible(NULL))
    if (!is.null(node$t) && node$t %in% c("Link", "Image")) {
      links <<- c(links, node$c[[3L]][[1L]])
    }
    if (!is.null(node$t) && node$t %in% c("RawInline", "RawBlock") && identical(node$c[[1L]], "html")) {
      attributes <- regmatches(node$c[[2L]], gregexpr(
        "\\b(href|src)\\s*=\\s*(\"[^\"]*\"|'[^']*'|[^[:space:]>]+)",
        node$c[[2L]], ignore.case = TRUE, perl = TRUE))[[1L]]
      targets <- sub("^[^=]+=\\s*", "", attributes, perl = TRUE)
      targets <- sub("^[\"'](.*)[\"']$", "\\1", targets)
      links <<- c(links, targets)
    }
    invisible(lapply(node, visit))
  }
  visit(document$blocks)
  links <- unique(links)
  parts <- tools::parse_URI_reference(links)
  local <- (parts$scheme == "file" | (parts$scheme == "" & !nzchar(parts$authority))) & nzchar(parts$path)
  links <- links[local]
  paths <- vapply(parts$path[local], utils::URLdecode, "")
  source <- normalizePath(source, winslash = "/", mustWork = TRUE)
  absolute <- parts$scheme[local] == "file" | grepl("^(/|[A-Za-z]:|\\\\)", paths)
  resolved <- vapply(file.path(source, paths), normalizePath, "",
                     winslash = "/", mustWork = FALSE)
  present <- !absolute & startsWith(resolved, paste0(source, "/")) & file.exists(resolved)
  report <- data.frame(URI = links, Path = paths, PresentInArchive = present)
  utils::write.csv(report, file.path(evidence, "archive-readme-links.csv"), row.names = FALSE)
  if (any(!present)) {
    stop("README links are absent from the exact source archive: ",
         paste(unique(paths[!present]), collapse = ", "))
  }
  invisible(report)
}

ftq_cran_unexpected_notes <- function(notes, maintainer) {
  # These are domain terms, not blanket exemptions for spelling diagnostics.
  spelling <- "^(ROI|Tractograms|tractograms|tractography) \\([0-9]+:[0-9]+(, [0-9]+:[0-9]+)*\\)$"
  allowed <- function(note) {
    lines <- trimws(strsplit(note, "\n", fixed = TRUE)[[1L]])
    lines <- chartr("\u2018\u2019", "''", lines[nzchar(lines)])
    if (!length(lines) || !grepl("^checking CRAN incoming feasibility \\.\\.\\.( \\[[^]]+\\])? NOTE$", lines[1L])) {
      return(FALSE)
    }
    body <- lines[-1L]
    heading <- body == "Possibly misspelled words in DESCRIPTION:"
    words <- grepl(spelling, body)
    if (any(heading) != any(words)) return(FALSE)
    any(body == "New submission" | words) &&
      all(body %in% c(paste0("Maintainer: '", maintainer, "'"), "New submission",
                      "Possibly misspelled words in DESCRIPTION:") | words)
  }
  notes[!vapply(notes, allowed, logical(1L))]
}

ftq_cran_environment <- function(profile, evidence, spelling = FALSE,
                                 subprocess = c("check", "direct")) {
  profile <- match.arg(profile, c("normal", "long-temp"))
  subprocess <- match.arg(subprocess)
  env <- c(NOT_CRAN = "false", `_R_CHECK_FORCE_SUGGESTS_` = "true",
           `_R_CHECK_CRAN_INCOMING_` = "true",
           `_R_CHECK_CRAN_INCOMING_REMOTE_` = "true",
           `_R_CHECK_CRAN_INCOMING_CHECK_FILE_URIS_` = "true",
           `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = if (spelling) "true" else "false",
           `_R_CHECK_RD_VALIDATE_RD2HTML_` = "true",
           `_R_CHECK_SYSTEM_CLOCK_` = "true",
           `_R_CHECK_TIMINGS_` = "0")
  windows <- .Platform$OS.type == "windows"
  if (profile == "long-temp" && !windows) stop("The long-temp profile requires Windows")
  if (windows && (subprocess == "check" || profile == "long-temp")) {
    root <- normalizePath(Sys.getenv("RUNNER_TEMP"), winslash = "/", mustWork = TRUE)
    # --as-cran adds /RtmpXXXXXX/working_dir/RtmpXXXXXX (34 characters).
    # A direct R child adds only /RtmpXXXXXX (11 characters).
    suffix <- if (subprocess == "check") 34L else 11L
    target <- if (profile == "long-temp") 46L else 44L
    parent_length <- max(nchar(root), target - suffix)
    if (parent_length == nchar(root) + 1L) parent_length <- parent_length + 1L
    if (parent_length + suffix > 50L) stop("Runner temp root is too long to calibrate CRAN-depth paths")
    padding <- parent_length - nchar(root)
    parent <- if (!padding) root else file.path(root, strrep("t", padding - 1L))
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    env <- c(env, TMPDIR = parent, TMP = parent, TEMP = parent)
  }
  # These variables must reach a new R process; changing the current tempdir is ineffective.
  probe <- callr::r(function() {
    path <- tempfile("write-probe-")
    writable <- file.create(path)
    unlink(path)
    list(StartupProbeTempdir = normalizePath(tempdir(), winslash = "/"),
         StartupProbeCharacters = nchar(normalizePath(tempdir(), winslash = "/")),
         Writable = writable, RVersion = R.version.string)
  }, env = env, user_profile = FALSE, system_profile = FALSE)
  expected <- probe$StartupProbeCharacters + if (subprocess == "check") 23L else 0L
  write.dcf(as.data.frame(c(list(Profile = profile, Subprocess = subprocess,
                                ExpectedTestTempdirCharacters = expected), probe, as.list(env)),
                          check.names = FALSE),
            file.path(evidence, "check-environment.dcf"))
  if (!isTRUE(probe$Writable)) stop("Fresh R process could not write its temporary directory")
  if (windows && (subprocess == "check" || profile == "long-temp") &&
      !(expected >= 44L && expected <= 50L)) {
    stop("Expected subprocess tempdir did not meet the requested 44-50 character depth")
  }
  env
}

ftq_cran_tempdir_hook <- function(evidence) {
  observations <- file.path(evidence, "process-tempdirs")
  dir.create(observations, recursive = TRUE, showWarnings = FALSE)
  observations <- normalizePath(observations, winslash = "/", mustWork = TRUE)
  hook <- file.path(evidence, "record-check-tempdir.R")
  expression <- substitute(local({
    process <- Sys.getpid()
    directory <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
    working <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    record <- data.frame(ProcessID = process, WorkingDirectory = working,
                          Tempdir = directory, Characters = nchar(directory),
                          Writable = file.access(directory, 2L) == 0L)
    write.dcf(record, file.path(OBSERVATIONS, paste0("process-", process, ".dcf")))
    if (grepl("[.]Rcheck/tests$", working)) {
      reg.finalizer(.GlobalEnv, function(...) {
        capture.output(sessionInfo(), file = file.path(OBSERVATIONS, paste0("session-", process, ".txt")))
      }, onexit = TRUE)
    }
  }), list(OBSERVATIONS = observations))
  writeLines(deparse(expression), hook)
  normalizePath(hook, winslash = "/", mustWork = TRUE)
}

ftq_cran_validate_test_tempdirs <- function(evidence, package, windows = .Platform$OS.type == "windows") {
  files <- list.files(file.path(evidence, "process-tempdirs"), "^process-.*[.]dcf$", full.names = TRUE)
  if (!length(files)) stop("No actual R CMD check subprocess tempdirs were recorded")
  records <- do.call(rbind, lapply(files, function(file) as.data.frame(read.dcf(file), stringsAsFactors = FALSE)))
  tests <- normalizePath(file.path(evidence, paste0(package, ".Rcheck"), "tests"),
                          winslash = "/", mustWork = TRUE)
  observed <- records[tolower(records$WorkingDirectory) == tolower(tests), , drop = FALSE]
  utils::write.csv(observed, file.path(evidence, "actual-test-tempdirs.csv"), row.names = FALSE)
  if (!nrow(observed)) stop("The actual package-test R process did not record its temporary directory")
  lengths <- as.integer(observed$Characters)
  if (anyNA(lengths) || !all(lengths == nchar(observed$Tempdir)) || !all(observed$Writable == "TRUE")) {
    stop("Actual test temporary-directory evidence is inconsistent or not writable")
  }
  if (windows && !all(lengths >= 44L & lengths <= 50L)) {
    stop("Actual R CMD check tests did not use the required 44-50 character temporary-directory depth")
  }
  invisible(observed)
}

ftq_cran_setup_aspell <- function(evidence) {
  dir.create(evidence, recursive = TRUE, showWarnings = FALSE)
  probe <- function(program) {
    if (!nzchar(program)) return(FALSE)
    output <- suppressWarnings(system2(program,
      c("--master=en_US", "--add-extra-dicts=en_GB", "list"),
      input = "colour color zzfibrtraqspellingprobezz", stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    (is.null(status) || identical(status, 0L)) && identical(output, "zzfibrtraqspellingprobezz")
  }
  aspell <- unname(Sys.which("aspell"))
  if (!probe(aspell)) {
    brew <- unname(Sys.which("brew"))
    if (!nzchar(brew)) brew <- "/home/linuxbrew/.linuxbrew/bin/brew"
    if (!file.exists(brew)) stop("Neither usable Aspell dictionaries nor runner Homebrew was found")
    Sys.setenv(HOMEBREW_NO_AUTO_UPDATE = "1")
    if (!identical(system2(brew, c("install", "aspell"), timeout = 600), 0L)) {
      stop("Runner-local Aspell installation failed")
    }
    prefix <- system2(brew, c("--prefix", "aspell"), stdout = TRUE)
    if (length(prefix) != 1L || !is.null(attr(prefix, "status"))) stop("Aspell prefix lookup failed")
    aspell <- file.path(prefix, "bin", "aspell")
  }
  if (!probe(aspell)) stop("Aspell could not verify both en_US and en_GB dictionaries")
  writeLines(c(system2(aspell, "--version", stdout = TRUE),
               system2(aspell, c("dump", "dicts"), stdout = TRUE)),
             file.path(evidence, "aspell-dictionaries.txt"))
  cat(dirname(aspell), "\n", sep = "", file = Sys.getenv("GITHUB_PATH"), append = TRUE)
  cat("FTQ_ASPELL_VERIFIED=true\n", file = Sys.getenv("GITHUB_ENV"), append = TRUE)
}

ftq_cran_report_failure <- function(result) {
  report_error <- length(result$error) == 1L &&
    grepl("Main report: (cannot open the connection|unable to start png\\(\\) device)", result$error)
  quoted <- unlist(regmatches(result$warnings, gregexpr("'[^']+'", result$warnings)))
  paths <- sub("'$", "", sub("^'", "", quoted))
  paths <- unique(paths[grepl("^[A-Za-z]:", paths) & nchar(paths) >= 260L])
  list(reproduced = report_error && length(paths) > 0L,
       recognized_report_error = report_error, long_paths = paths,
       path_characters = nchar(paths))
}

ftq_cran_baseline_disposition <- function(result, configuration) {
  diagnosis <- ftq_cran_report_failure(result)
  files <- result$files
  measured <- is.data.frame(files) && nrow(files) > 0L &&
    all(c("Source", "Characters", "Copied") %in% names(files)) &&
    is.character(files$Source) && is.numeric(files$Characters) &&
    is.logical(files$Copied) && all(!is.na(files$Source)) &&
    isTRUE(all(files$Characters == nchar(files$Source))) && isTRUE(all(files$Copied))
  maximum <- if (measured) max(files$Characters) else NA_integer_
  enabled <- is.numeric(configuration$LongPathsEnabled) &&
    length(configuration$LongPathsEnabled) == 1L &&
    isTRUE(configuration$LongPathsEnabled == 1L)
  completed <- is.null(result$error) && !length(result$warnings) &&
    isTRUE(result$example_passed)
  outcome <- if (diagnosis$reproduced) {
    "reproduced"
  } else if (completed && enabled && isTRUE(maximum >= 260L)) {
    "not_reproduced_long_paths_enabled"
  } else {
    "unexpected_outcome"
  }
  c(diagnosis, list(outcome = outcome, accepted = outcome != "unexpected_outcome",
                   long_paths_enabled = enabled, maximum_generated_path = maximum))
}
