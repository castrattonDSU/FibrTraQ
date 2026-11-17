ftq_preflight_manifest <- function(m) {
  paths <- c(
    m$fnode, m[["local"]], m$tractogram, m$target_paths,
    m$endpoint_context_paths, m$ref_image, ftq_expected_paths(m$expected)
  )
  missing <- paths[!is.na(paths) & !file.exists(paths)]
  if (length(missing)) stop("Input file(s) do not exist:\n", paste(missing, collapse = "\n"))
  invisible(TRUE)
}
