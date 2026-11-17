ftq_unified_example_summary <- function(validation) {
  if (is.null(validation)) return("")
  esc <- ftq_report_html_escape
  count <- function(x) paste0(sum(x$passed, na.rm = TRUE), "/", nrow(x))
  details <- vapply(c("tables", "claims", "reports"), function(family) {
    checks <- validation[[family]]
    label <- if (family == "claims") checks$claim else checks$artifact
    items <- paste0("<li>", esc(label), ": ",
                    ifelse(checks$passed, "passed", "FAILED"),
                    " &mdash; ", esc(checks$details), "</li>")
    paste0("<h3>", family, "</h3><ul>", paste(items, collapse = ""), "</ul>")
  }, character(1))
  paste0("<details id='ftq-example-checks'><summary>Bundled example verification</summary>",
         "<p>Known-answer comparisons: ", count(validation$tables), " tables and ",
         count(validation$claims), " claims passed. Computational validation: ",
         if (isTRUE(validation$package_validation_passed)) "passed" else "not passed",
         ".</p>", paste(details, collapse = "\n"), "</details>")
}

ftq_unified_example_validation <- function(path, validation) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  at <- match("<!-- FTQ_EXAMPLE_VERIFICATION -->", lines)
  if (is.na(at)) stop("Cannot locate the example verification section in ", path)
  lines[at] <- ftq_unified_example_summary(validation)
  ftq_write_report_page(lines, path)
}
