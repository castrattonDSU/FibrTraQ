ftq_write_regional_study_index <- function(run_index,
                                           path,
                                           output_dir,
                                           description,
                                           example_validation = NULL) {
  esc <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "unavailable"
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    gsub("'", "&#39;", x, fixed = TRUE)
  }
  field <- function(row, name, default = "unavailable") {
    x <- row[[name]]
    if (is.null(x) || !length(x) || is.na(x[1L]) || !nzchar(as.character(x[1L]))) {
      return(default)
    }
    as.character(x[1L])
  }
  link <- function(target, label) {
    if (!length(target) || is.na(target) || !nzchar(target)) return("")
    href <- ftq_report_href(target, dirname(ftq_report_normalize_path(path)))
    paste0("<a href='", esc(href), "'>", esc(label), "</a>")
  }
  mode_text <- function(row) {
    mode <- field(row, "target_intersection_mode", "not recorded")
    labels <- c(full_streamline = "Full streamline", endpoint_voxels = "Endpoint voxels",
                local_search = "Local search")
    text <- if (mode %in% names(labels)) paste0(labels[[mode]], " (", mode, ")") else mode
    if (identical(mode, "local_search")) {
      radius <- suppressWarnings(as.numeric(field(row, "local_search_radius_mm", NA_character_)))
      text <- paste0(text, "; radius ", if (is.finite(radius))
        paste0(format(radius, trim = TRUE), " mm") else "unavailable")
    }
    allocation <- field(row, "focal_assignment_mode", "")
    if (nzchar(allocation)) text <- paste0(text, "; focal allocation: ", allocation)
    text
  }
  modes <- vapply(seq_len(nrow(run_index)), function(i) mode_text(run_index[i, , drop = FALSE]),
                  character(1))
  shared_mode <- length(unique(modes)) == 1L
  download_labels <- c(
    result_csv = "Complete per-voxel data (CSV)",
    voxel_data_csv = "Complete per-voxel data with coordinates (CSV)",
    rsc_csv = "Complete RSC summary (CSV)",
    fate_csv = "Voxel-weighted fate (CSV)",
    roi_summary_csv = "Detailed run summary and unique-streamline fate (CSV)",
    endpoint_summary_csv = "Endpoint assignment audit (CSV)",
    endpoint_penetration_csv = "Target penetration (CSV)",
    endpoint_context_csv = "Endpoint coverage (CSV)",
    local_search_assignments_csv = "Local-search assignments (CSV)",
    metadata_json = "Methods and provenance (JSON)",
    validation_json = "Computational validation (JSON)"
  )
  rows <- vapply(seq_len(nrow(run_index)), function(i) {
    row <- run_index[i, , drop = FALSE]
    legacy_completed <- identical(field(row, "run_status"), "completed")
    analysis <- field(row, "analysis_status", if (legacy_completed) "completed" else "unavailable")
    validation <- field(row, "validation_status", if (isTRUE(row$validation_passed)) "passed" else "not passed")
    primary <- field(row, "primary_report", "")
    action <- if (nzchar(primary)) link(primary, "Results") else paste0(
      "<a href='#run-", i, "'>Results and available files</a>"
    )
    auxiliary <- c(link(field(row, "exploratory_report", ""), "Target report"),
                   link(field(row, "interactive_3d_qc", ""), "3D QC"))
    auxiliary <- auxiliary[nzchar(auxiliary)]
    downloads <- vapply(names(download_labels), function(name) {
      target <- field(row, name, "")
      if (!nzchar(target)) return("")
      # Historical summary CSVs can contain original absolute paths after a
      # move. Do not advertise stale download links; repaired reports remain.
      if (!file.exists(target)) return("")
      paste0("<li>", link(target, download_labels[[name]]), "</li>")
    }, character(1))
    status_fields <- c(primary_report = "Main report", exploratory_report = "Target report",
                       interactive_3d_qc = "3D QC")
    presentation <- vapply(names(status_fields), function(name) {
      status <- field(row, paste0(name, "_status"), "")
      if (!nzchar(status)) return("")
      paste0(esc(status_fields[[name]]), ": ", esc(gsub("_", " ", status, fixed = TRUE)))
    }, character(1))
    presentation <- presentation[nzchar(presentation)]
    error <- field(row, "error", "")
    detail <- paste0(
      "<details id='run-", i, "'", if (nzchar(error) || !nzchar(primary)) " open" else "", ">",
      "<summary>Supporting views, downloads and status</summary>",
      if (length(presentation)) paste0("<p>", paste(presentation, collapse = "; "), ".</p>") else "",
      if (length(auxiliary)) paste0("<p>", paste(auxiliary, collapse = " &middot; "), "</p>") else "",
      if (any(nzchar(downloads))) paste0("<ul>", paste(downloads, collapse = ""), "</ul>") else
        "<p class='muted'>No additional data links are available in this index.</p>",
      "</details>"
    )
    paste0(
      "<tr><th scope='row'>", esc(row$subject), "</th><td>", esc(row$focal_region), "</td>",
      "<td class='", if (analysis == "completed") "pass" else "fail", "'>", esc(analysis), "</td>",
      "<td class='", if (validation == "passed") "pass" else "fail", "'>", esc(validation), "</td>",
      "<td>", esc(gsub("_", " ", field(row, "biological_qc_status"), fixed = TRUE)), "</td>",
      if (!shared_mode) paste0("<td>", esc(modes[[i]]), "</td>") else "",
      "<td>", action, detail, "</td></tr>",
      if (nzchar(error)) paste0("<tr><td colspan='", if (shared_mode) 6L else 7L,
        "' class='fail'>Attention: ", esc(error), "</td></tr>") else ""
    )
  }, character(1))
  title <- if (nrow(run_index) == 1L) {
    paste0("FibrTraQ Results: ", run_index$subject[[1L]], " / ", run_index$focal_region[[1L]])
  } else "FibrTraQ Results"
  guide <- file.path(output_dir, "OUTPUTS_README.html")
  root_files <- c("fibrtraq-run-index.csv" = "Run index (CSV)",
                  "fibrtraq-prepared-inputs.csv" = "Prepared-input provenance (CSV)")
  root_links <- vapply(names(root_files), function(name) {
    target <- file.path(output_dir, name)
    if (file.exists(target)) link(target, root_files[[name]]) else ""
  }, character(1))
  root_links <- c(root_links[nzchar(root_links)],
                  if (file.exists(guide)) link(guide, "Output definitions and guide"))
  example <- character()
  if (!is.null(example_validation)) {
    counts <- function(x) paste0(sum(x$passed, na.rm = TRUE), "/", nrow(x))
    checked <- isTRUE(example_validation$package_validation_passed) &&
      nrow(example_validation$tables) > 0L && nrow(example_validation$claims) > 0L &&
      isTRUE(all(example_validation$tables$passed)) && isTRUE(all(example_validation$claims$passed))
    example <- c(
      "<section><h2>Bundled example verification</h2>",
      paste0("<p class='", if (checked) "pass" else "fail", "'>Known-answer comparisons: ",
        counts(example_validation$tables), " tables and ", counts(example_validation$claims),
        " claims passed. Computational validation: ",
        if (isTRUE(example_validation$package_validation_passed)) "passed" else "not passed", ".</p>"),
      paste0("<p>Generated HTML files present and nonempty: ", counts(example_validation$reports),
        ". This checks file generation; browser rendering and interactive behavior require separate review.</p>"),
      "<details><summary>Example comparison details</summary>"
    )
    for (family in c("tables", "claims", "reports")) {
      checks <- example_validation[[family]]
      label <- if (family == "claims") checks$claim else checks$artifact
      example <- c(example, paste0("<h3>", esc(family), "</h3><ul>"),
        vapply(seq_len(nrow(checks)), function(i) paste0("<li>", esc(label[[i]]), ": ",
          if (isTRUE(checks$passed[[i]])) "passed" else "FAILED", " \u2014 ",
          esc(checks$details[[i]]), "</li>"), character(1)), "</ul>")
    }
    example <- c(example, "</details></section>")
  }
  html <- c(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    paste0("<title>", esc(title), "</title>"),
    paste0("<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;",
      "max-width:1200px;margin:32px auto;padding:0 20px;color:#1f2933;line-height:1.5;overflow-wrap:anywhere}",
      ".table-scroll{width:100%;overflow-x:auto;-webkit-overflow-scrolling:touch}",
      ".table-scroll:focus{outline:2px solid #2563eb;outline-offset:2px}",
      "table{border-collapse:collapse;width:100%;min-width:700px}",
      "th,td{border:1px solid #d8dee4;padding:10px;text-align:left;vertical-align:top}",
      "thead th{background:#f3f4f6}tbody th{font-weight:600}.pass{color:#176b35}",
      ".fail{color:#9f1239;font-weight:700}.muted{color:#596579}a{color:#175db0}",
      "details{margin-top:10px}summary{cursor:pointer}ul{padding-left:20px}",
      "@media(max-width:700px){body{margin:16px 0;padding:0 12px}h1{font-size:1.6rem}",
      "th,td{padding:6px 8px;font-size:.9rem}}</style></head><body>"),
    paste0("<h1>", esc(title), "</h1>"),
    paste0("<p>", nrow(run_index), if (nrow(run_index) == 1L) " run recorded. " else " runs recorded. ",
      sum(run_index$validation_passed, na.rm = TRUE), " passed computational validation.</p>"),
    if (!nrow(run_index)) "<p>Inputs are being prepared; no analysis run has finished yet.</p>" else "",
    if (shared_mode) paste0("<p><strong>Method:</strong> ", esc(modes[[1L]]), ".</p>") else "",
    "<p>Computational validation checks numerical consistency. Biological QC is a separate review of anatomical context and alignment; unavailable or unreviewed QC is not evidence of anatomical validity.</p>",
    example,
    "<div class='table-scroll' role='region' tabindex='0'>",
    "<table><thead><tr><th>Subject</th><th>Focal region</th><th>Analysis</th><th>Validation</th><th>Biological QC</th>",
    if (!shared_mode) "<th>Method</th>" else "",
    "<th>Results</th></tr></thead><tbody>", rows, "</tbody></table></div>",
    if (length(root_links)) paste0("<h2>Files and definitions</h2><p>",
      paste(root_links, collapse = " &middot; "), "</p>") else "",
    "<p class='muted'>Preparation records retain optional inputs as unavailable when none were supplied. A missing optional reference or b-value file is not an analysis failure. Keep the complete results folder together and extract it before opening the HTML.</p>",
    "</body></html>"
  )
  writeLines(enc2utf8(html), path, useBytes = TRUE)
  invisible(path)
}

ftq_regional_study_relative_path <- function(path, root) {
  ftq_report_relative_path(path, root)
}
