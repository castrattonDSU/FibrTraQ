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

ftq_cran_environment <- function(profile, evidence, spelling = FALSE) {
  profile <- match.arg(profile, c("normal", "long-temp"))
  env <- c(NOT_CRAN = "false", `_R_CHECK_FORCE_SUGGESTS_` = "true",
           `_R_CHECK_CRAN_INCOMING_` = "true",
           `_R_CHECK_CRAN_INCOMING_REMOTE_` = "true",
           `_R_CHECK_CRAN_INCOMING_CHECK_FILE_URIS_` = "true",
           `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = if (spelling) "true" else "false",
           `_R_CHECK_RD_VALIDATE_RD2HTML_` = "true",
           `_R_CHECK_SYSTEM_CLOCK_` = "true")
  if (profile == "long-temp") {
    if (.Platform$OS.type != "windows") stop("The long-temp profile requires Windows")
    root <- normalizePath(Sys.getenv("RUNNER_TEMP"), winslash = "/", mustWork = TRUE)
    # R appends /RtmpXXXXXX: target 46 characters (CRAN's observed prefix was 44).
    padding <- 46L - 11L - nchar(root) - 1L
    if (padding < 4L) stop("Runner temp root is too long to calibrate CRAN-depth paths")
    parent <- file.path(root, paste0("cran", strrep("t", padding - 4L)))
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    env <- c(env, TMPDIR = parent, TMP = parent, TEMP = parent)
  }
  # These variables must reach a new R process; changing the current tempdir is ineffective.
  probe <- callr::r(function() {
    path <- tempfile("write-probe-")
    writable <- file.create(path)
    unlink(path)
    list(Tempdir = normalizePath(tempdir(), winslash = "/"),
         TempdirCharacters = nchar(normalizePath(tempdir(), winslash = "/")),
         Writable = writable, RVersion = R.version.string)
  }, env = env, user_profile = FALSE, system_profile = FALSE)
  write.dcf(as.data.frame(c(list(Profile = profile), probe, as.list(env)), check.names = FALSE),
            file.path(evidence, "check-environment.dcf"))
  if (!isTRUE(probe$Writable)) stop("Fresh R process could not write its temporary directory")
  if (profile == "long-temp" && !(probe$TempdirCharacters >= 44L && probe$TempdirCharacters <= 50L)) {
    stop("Fresh R tempdir did not reproduce the requested 44-50 character depth")
  }
  env
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
