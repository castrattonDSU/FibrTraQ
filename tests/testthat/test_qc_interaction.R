test_that("generated QC JavaScript supports stable trackball interaction", {
  node <- unname(Sys.which("node"))
  skip_if(!nzchar(node), "Node.js is unavailable for the generated viewer interaction test")
  harness <- test_path("..", "js", "qc-interaction.cjs")
  expect_true(file.exists(harness))
  html <- tempfile("qc-interaction-", fileext = ".html")
  on.exit(unlink(html), add = TRUE)

  # Synthetic scene coordinates deliberately have a nonzero centre so the
  # JavaScript checks exercise both the RAS basis and centring transformation.
  points <- list(c(12, -7, 5), c(14, -6, 8), c(10, -8, 3))
  scene <- list(
    metadata = list(title = "Synthetic interaction fixture", subtitle = "Test only",
      interaction = "Drag to rotate; scroll or pinch to zoom.", method = "Fixture",
      context = "Synthetic coordinates", note = "No biological data"),
    bounds = list(center = c(12, -7, 5), radius = 10),
    brain = points, focal = points,
    targets = list(list(target = "Synthetic destination", color = "#1479b8",
      point = c(14, -6, 8), boundary = points)),
    legend = list(list(target = "Synthetic destination", color = "#1479b8")),
    endpointAnalysis = FALSE, endpointClasses = list(),
    streamlines = list(list(color = "#1479b8", points = points, focal = points,
      targetSegment = points, after = points, hits = list(points[[1L]]),
      targetHits = list(points[[2L]])))
  )
  writeLines(FibrTraQ:::ftq_qc_html_template(
    jsonlite::toJSON(scene, auto_unbox = TRUE, digits = 16)
  ), html)
  output <- suppressWarnings(system2(node,
    c(shQuote(normalizePath(harness, winslash = "/", mustWork = TRUE)),
      shQuote(normalizePath(html, winslash = "/", mustWork = TRUE))),
    stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_identical(as.integer(status), 0L, info = paste(output, collapse = "\n"))
  expect_true(any(grepl("All generated QC interaction checks passed", output,
                       fixed = TRUE)), info = paste(output, collapse = "\n"))
})
