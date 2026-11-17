# In a checkout matching your installed package:
#   research <- source("archive/load.R")$value
local({
  source_calls <- which(vapply(seq_along(sys.frames()), function(i) {
    identical(sys.function(i), base::source)
  }, logical(1)))
  if (!length(source_calls)) stop("Load this file with source().", call. = FALSE)
  caller <- sys.frame(tail(source_calls, 1L))
  directory <- if (isTRUE(caller$chdir)) getwd() else dirname(caller$ofile)

  core <- asNamespace("FibrTraQ")
  if (!all(c("fibrtraq_single", "fibrtraq_study") %in% getNamespaceExports(core))) {
    stop("Install FibrTraQ from the same revision as this archive.", call. = FALSE)
  }
  research <- new.env(parent = core)
  research$.archive_root <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  for (path in list.files(file.path(research$.archive_root, "R"), pattern = "[.]R$", full.names = TRUE)) {
    sys.source(path, envir = research)
  }

  # These former public commands are still called by the installed workflows.
  # Expose their actual implementations here for historical research scripts.
  core_helpers <- c(
    "FibrTraQ", "as_track", "fibrtraq", "ftq_endpoint_tables",
    "ftq_exploratory_report", "ftq_exploratory_summary", "ftq_interactive_qc_viewer",
    "ftq_local_search_tables", "ftq_make_union_mask", "ftq_preflight", "ftq_read_manifest",
    "ftq_report", "ftq_run", "ftq_validate_run", "ftq_write_outputs",
    "iso_coords", "rsc_summary", "summary_connectivity"
  )
  list2env(mget(core_helpers, envir = core, inherits = FALSE), envir = research)
  research$ftq_run_files <- core$fibrtraq_single
  research$ftq_run_regional_study <- core$fibrtraq_study
  # Scripts using the archived 0.3.5 study command used a regional union by default.
  formals(research$ftq_run_regional_study)$local_masks <- "union"
  registerS3method("print", "ftq_example_result", research$print.ftq_example_result,
                   envir = asNamespace("base"))
  research
})
