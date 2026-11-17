# Only the explicit navigation region in a freshly generated primary report is
# eligible for replacement. Historical HTML without this marker is untouched.
ftq_report_add_navigation <- function(html, links) {
  if (!is.character(html) || length(html) != 1L || is.na(html) ||
      !file.exists(html) || dir.exists(html)) {
    stop("html must identify an existing report file.", call. = FALSE)
  }
  text <- rawToChar(readBin(html, "raw", n = file.info(html)$size))
  Encoding(text) <- "UTF-8"
  start_marker <- "<!-- FIBRTRAQ_NAVIGATION -->"
  end_marker <- "<!-- /FIBRTRAQ_NAVIGATION -->"
  starts <- gregexpr(start_marker, text, fixed = TRUE)[[1L]]
  ends <- gregexpr(end_marker, text, fixed = TRUE)[[1L]]
  if (starts[[1L]] == -1L && ends[[1L]] == -1L) return(invisible(FALSE))
  if (length(starts) != 1L || length(ends) != 1L ||
      starts[[1L]] < 0L || ends[[1L]] < starts[[1L]]) {
    stop("The report navigation marker is missing, repeated or out of order.",
         call. = FALSE)
  }
  if (length(links) && (is.null(names(links)) || any(!nzchar(names(links))))) {
    stop("Navigation links must have readable names.", call. = FALSE)
  }
  root <- dirname(html)
  anchors <- vapply(seq_along(links), function(i) {
    path <- links[[i]]
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !file.exists(path) || dir.exists(path)) return("")
    paste0("<a href='", ftq_report_html_escape(ftq_report_href(path, root)),
           "'>", ftq_report_html_escape(names(links)[[i]]), "</a>")
  }, character(1))
  anchors <- anchors[nzchar(anchors)]
  content <- if (length(anchors)) {
    paste0("<nav class='report-navigation' aria-label='Related reports'>",
           paste(anchors, collapse = " &middot; "), "</nav>")
  } else ""
  replacement <- paste0(start_marker, content, end_marker)
  updated <- paste0(
    substr(text, 1L, starts[[1L]] - 1L), replacement,
    substring(text, ends[[1L]] + nchar(end_marker))
  )
  if (!identical(text, updated)) writeBin(charToRaw(enc2utf8(updated)), html)
  invisible(TRUE)
}
