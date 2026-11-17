# The original renderers produce escaped HTML with links relative to the final
# report. Reuse their scientific tables and figures in a single document.
ftq_report_section_body <- function(document) {
  if (!length(document)) return("")
  html <- paste(document, collapse = "\n")
  if (!nzchar(html)) return("")
  start <- regexpr("<body>", html, fixed = TRUE)[1L]
  finish <- regexpr("</body>", html, fixed = TRUE)[1L]
  if (start < 0L || finish <= start) stop("Incomplete report section.")
  body <- substr(html, start + nchar("<body>"), finish - 1L)
  body <- sub("^\\s*<h1>[^<]*</h1>", "", body, perl = TRUE)
  gsub("<img ", "<img loading='lazy' ", body, fixed = TRUE)
}

ftq_write_report_page <- function(lines, path) {
  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  staged <- tempfile(".ftq-", tmpdir = directory)
  on.exit(unlink(staged), add = TRUE)
  writeLines(enc2utf8(lines), staged, useBytes = TRUE)
  if (!file.rename(staged, path)) {
    stop("Could not replace the report at ", path, ". Check the path and write access.")
  }
  invisible(path)
}

ftq_write_unified_report <- function(path, output_dir, run_index, sections,
                                     description, planned_runs = nrow(run_index)) {
  esc <- ftq_report_html_escape
  root <- dirname(path)
  n <- nrow(run_index)
  labels <- if (n) paste(run_index$subject, run_index$focal_region, sep = " / ") else character()
  failed <- run_index$run_status %in% "failed"
  failed_count <- sum(failed)
  incomplete <- n < planned_runs
  study_notice <- if (incomplete || (planned_runs > 1L && failed_count > 0L)) {
    paste0("<p id='ftq-study-status' class='warning' role='status'><strong>Study incomplete.</strong>",
      if (incomplete) paste0(" Run records available for ", n, " of ", planned_runs, " planned analyses."),
      if (failed_count == 1L) " One analysis failed.",
      if (failed_count > 1L) paste0(" ", failed_count, " analyses failed."),
      if (failed_count && n > 1L) " Select a run marked (failed) for details.",
      "</p>")
  }
  root_paths <- file.path(output_dir, c("fibrtraq-run-index.csv", "fibrtraq-prepared-inputs.csv"))
  downloads <- ftq_report_downloads(as.list(root_paths), root,
                                    c("Run index (CSV)", "Prepared inputs (CSV)"))
  field <- function(row, name) {
    value <- row[[name]]
    if (!length(value) || is.na(value[1L])) "" else as.character(value[1L])
  }
  run_html <- function(i) {
    row <- run_index[i, , drop = FALSE]
    entry <- sections[[i]]
    primary <- ftq_report_section_body(entry$primary)
    target <- ftq_report_section_body(entry$target)
    scene <- entry$viewer
    has_scene <- length(scene) == 1L && !is.na(scene) && file.exists(scene)
    stages <- c(analysis_status = "Analysis", validation_status = "Computational validation",
                biological_qc_status = "Biological QC", primary_report_status = "Results",
                exploratory_report_status = "Target figures", interactive_3d_qc_status = "3D viewer")
    statuses <- vapply(names(stages), function(name) {
      value <- field(row, name)
      if (!nzchar(value)) return("")
      paste0("<li>", stages[[name]], ": ", esc(gsub("_", " ", value, fixed = TRUE)), "</li>")
    }, character(1))
    error <- field(row, "error")
    run_files <- c(result_csv = "Per-voxel data", rsc_csv = "RSC summary", fate_csv = "Fate summary",
                   roi_summary_csv = "ROI summary", endpoint_summary_csv = "Endpoint summary",
                   endpoint_context_csv = "Endpoint context", endpoint_penetration_csv = "Endpoint penetration",
                   local_search_assignments_csv = "Local-search assignments",
                   metadata_json = "Run metadata (JSON)", validation_json = "Computational validation (JSON)")
    paths <- lapply(names(run_files), function(name) field(row, name))
    files <- ftq_report_downloads(paths, root, unname(run_files))
    paste0(
      "<h1 id='ftq-current-run'>", esc(labels[i]), "</h1>",
      if (nzchar(error)) paste0("<p class='warning' role='alert'>", esc(error), "</p>"),
      "<nav aria-label='Report sections'>",
      if (nzchar(primary)) "<a href='#ftq-primary'>Results</a>",
      if (nzchar(target)) "<a href='#ftq-targets'>Target figures</a>",
      if (has_scene) "<a href='#ftq-qc'>3D QC viewer</a>",
      "<a href='#ftq-files'>Run details and files</a></nav>",
      if (nzchar(primary)) paste0("<section id='ftq-primary' aria-label='Results'>", primary, "</section>"),
      if (nzchar(target)) paste0("<section id='ftq-targets'><h2>Target figures</h2>", target, "</section>"),
      if (has_scene) paste0("<section id='ftq-qc'><h2>Interactive 3D QC</h2>",
        "<button id='ftq-open-viewer' type='button' aria-controls='ftq-viewer-host' data-scene='",
        esc(ftq_report_href(scene, root)), "'>Open 3D viewer</button>",
        "<div id='ftq-viewer-host' aria-live='polite'></div></section>"),
      "<section id='ftq-files'><h2>Run details and files</h2><ul>",
      paste(statuses, collapse = ""), "</ul>", files, "</section>"
    )
  }
  option_labels <- labels
  option_labels[failed] <- paste(labels[failed], "(failed)")
  options <- if (n) paste0("<option value='", seq_len(n) - 1L, "'>", esc(option_labels), "</option>")
  selector <- if (n) paste0("<label", if (n == 1L) " hidden",
    ">Subject / focal region <select id='ftq-run'>", paste(options, collapse = ""), "</select></label>")
  template <- function(i) {
    paste0("<template id='ftq-run-", i - 1L, "'>",
           if (i > 1L) run_html(i) else "", "</template>")
  }
  includes_viewer <- any(vapply(sections, function(s) length(s$viewer) == 1L &&
      !is.na(s$viewer) && file.exists(s$viewer), logical(1)))
  shell <- if (includes_viewer) paste0(
    "<template id='ftq-viewer-shell'>",
    esc(ftq_qc_html_template("FTQ_EMBEDDED_SCENE_JSON")), "</template>"
  )
  script <- system.file("report", "unified.js", package = "FibrTraQ", mustWork = TRUE)
  style <- system.file("report", "unified.css", package = "FibrTraQ", mustWork = TRUE)
  lines <- c(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    '<meta name="fibrtraq-report-layout" content="unified-v1">',
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    paste0("<title>FibrTraQ results", if (n == 1L) paste0(" | ", esc(labels)), "</title>"),
    "<style>", readLines(style, warn = FALSE), "</style></head><body>",
    "<header><strong>FibrTraQ results</strong>", selector, "</header>", study_notice,
    "<main id='ftq-active'>",
    if (n) run_html(1L) else "<p>No completed analyses are available yet.</p>",
    "</main>", vapply(seq_len(n), template, character(1)), shell,
    "<footer><p>To share this report, zip the entire output folder containing ",
    "<code>fibrtraq-run-index.csv</code>. Extract the folder before opening the HTML.</p>",
    downloads, "</footer>",
    "<!-- FTQ_EXAMPLE_VERIFICATION -->",
    "<script>", readLines(script, encoding = "UTF-8", warn = FALSE), "</script>",
    "</body></html>"
  )
  ftq_write_report_page(lines, path)
}
