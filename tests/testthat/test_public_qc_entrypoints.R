ftq_public_scene <- function(path) {
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  block <- regmatches(html, regexpr(
    '<script id="scene-data" type="application/json">.*?</script>', html
  ))
  json <- sub('^<script id="scene-data" type="application/json">', '', block)
  list(html = html, payload = jsonlite::fromJSON(
    sub('</script>$', '', json), simplifyVector = FALSE
  ))
}

test_that("the public example distinguishes exact evidence from HTML presence", {
  example <- ftq_run_example(backend = "r", verbose = FALSE)
  expect_true(example$validation$passed)
  expect_identical(nrow(example$validation$tables), 7L)
  expect_identical(nrow(example$validation$claims), 18L)
  expect_identical(nrow(example$validation$reports), 4L)
  expect_true(all(example$validation$tables$passed))
  expect_true(all(example$validation$claims$passed))
  expect_true(all(example$validation$reports$passed))
  console <- paste(capture.output(print(example)), collapse = "\n")
  expect_match(console, "Exact tables: 7/7", fixed = TRUE)
  expect_match(console, "Known answers: 18/18", fixed = TRUE)
  expect_match(console, "HTML files:   4/4 generated and non-empty", fixed = TRUE)
  expect_match(console, "browser rendering not checked", fixed = TRUE)
  expect_match(console, "Biological validity is not assessed", fixed = TRUE)
  primary <- paste(readLines(example$primary_report, warn = FALSE), collapse = "\n")
  navigation <- regmatches(primary, regexpr("<!-- FIBRTRAQ_NAVIGATION -->.*<!-- /FIBRTRAQ_NAVIGATION -->", primary))
  expect_match(navigation, "Results index", fixed = TRUE)
  expect_match(navigation, "Target report", fixed = TRUE)
  expect_match(navigation, "3-D QC viewer", fixed = TRUE)
  index <- paste(readLines(example$report_index, warn = FALSE), collapse = "\n")
  expect_match(index, "Known-answer comparisons", fixed = TRUE)
  expect_match(index, "18", fixed = TRUE)
  expect_match(index, "non.?empty")

  expect_identical(basename(dirname(example$interactive_3d_qc)), "report-figures")
  expect_false(dir.exists(file.path(dirname(example$primary_report), "figures")))
  expect_true(file.exists(file.path(dirname(example$interactive_3d_qc), "OUTPUTS_README.html")))
  scene <- ftq_public_scene(example$interactive_3d_qc)
  # The fixture includes the diagonal [0, 0, 3] to [8, 8, 3]. A radius of 4
  # (half its longest box side) crops it after rotation in a narrow canvas.
  expect_equal(unlist(scene$payload$bounds$center), c(4, 4, 3))
  expect_equal(scene$payload$bounds$radius, sqrt(32), tolerance = 1e-6)
  scene_points <- c(scene$payload$brain, scene$payload$focal)
  for (target in scene$payload$targets) {
    scene_points <- c(scene_points, list(target$point), target$boundary)
  }
  for (streamline in scene$payload$streamlines) {
    for (mode in c("points", "focal", "targetSegment", "after", "hits", "targetHits")) {
      scene_points <- c(scene_points, streamline[[mode]])
    }
    scene_points <- c(scene_points, lapply(streamline$endpoints, `[[`, "point"))
  }
  distances <- vapply(scene_points, function(point) {
    sqrt(sum((unlist(point) - unlist(scene$payload$bounds$center))^2))
  }, numeric(1))
  expect_true(all(distances <= scene$payload$bounds$radius + 1e-6))
  expect_identical(scene$payload$samplingStrategy, "target_quota_scan_order")
  expect_identical(scene$payload$targetIntersectionMode, "full_streamline")
  expect_match(scene$payload$metadata$subtitle,
               "not a random or prevalence-weighted sample", fixed = TRUE)
  expect_match(scene$payload$metadata$method, "full_streamline", fixed = TRUE)
  expect_match(scene$payload$metadata$note,
               "continuous streamline endpoints", fixed = TRUE)
  expect_match(scene$payload$metadata$note,
               "After Focal Hit follows stored path order", fixed = TRUE)
  expect_lt(nchar(scene$payload$metadata$interaction), 160L)
  expect_match(scene$html, 'elHud.textContent=data.metadata.interaction', fixed = TRUE)
  expect_match(scene$html, '<details class="viewer-details" id="viewerDetails">', fixed = TRUE)
  expect_false(grepl('<details class="viewer-details" id="viewerDetails" open', scene$html, fixed = TRUE))
  expect_match(scene$html, 'id="endpointStatus" aria-live="polite"', fixed = TRUE)
  expect_match(scene$html, "circle = first stored point; square = last stored point", fixed = TRUE)
  expect_match(scene$html, "Stored order does not show biological direction", fixed = TRUE)
  expect_match(scene$html, "updateEndpointStatus();render()", fixed = TRUE)
  expect_true(all(vapply(scene$payload$streamlines, function(s) {
    all(vapply(s$endpoints, function(e) {
      e$end %in% c("start", "end", "single") &&
        e$statusExact %in% c("focal", "target", "shared", "grid_outside", "context_outside", "other") &&
        length(e$point) == 3L
    }, logical(1)))
  }, logical(1))))
})

test_that("the dataset opener accepts an index without historical companions", {
  root <- tempfile("public-open-index-")
  output_dir <- file.path(root, "outputs", "harvard-oxford-left-thalamus")
  dir.create(output_dir, recursive = TRUE)
  index <- file.path(output_dir, "fibrtraq-report-index.html")
  writeLines("<!doctype html><title>Existing results</title>", index)
  original <- unname(tools::md5sum(index))
  original_root <- Sys.getenv("FIBRTRAQ_DTOR_ROOT", unset = NA_character_)
  on.exit(if (is.na(original_root)) Sys.unsetenv("FIBRTRAQ_DTOR_ROOT") else
            Sys.setenv(FIBRTRAQ_DTOR_ROOT = original_root), add = TRUE)
  Sys.setenv(FIBRTRAQ_DTOR_ROOT = root)
  opened <- character()
  local_mocked_bindings(
    ftq_qc_open_output = function(path) { opened <<- c(opened, path); invisible(path) },
    .package = "FibrTraQ"
  )
  script <- system.file("scripts", "fibrtraq-dtor985-view-full-report.R", package = "FibrTraQ")
  expect_true(nzchar(script))
  env <- new.env(parent = globalenv())
  expect_output(source(script, local = env), "Supporting tables are incomplete")
  expect_identical(opened, index)
  expect_identical(unname(tools::md5sum(index)), original)
  expect_identical(list.files(output_dir), "fibrtraq-report-index.html")

  portable <- file.path(output_dir, "fibrtraq-report-index-portable.html")
  writeLines("<!doctype html><title>Portable results</title>", portable)
  opened <- character()
  invisible(capture.output(source(script, local = new.env(parent = globalenv()))))
  expect_identical(opened, portable)
})

test_that("the dataset recipe resolves an installed companion without a checkout", {
  script <- system.file("scripts", "fibrtraq-dtor985-start-to-finish.R", package = "FibrTraQ")
  expressions <- as.list(parse(script))
  assignment <- Filter(function(x) {
    is.call(x) && identical(x[[1L]], as.name("<-")) &&
      identical(x[[2L]], as.name("streamline_figure_script"))
  }, expressions)
  expect_length(assignment, 1L)
  env <- new.env(parent = environment())
  env$repo_root <- tempfile("absent-source-checkout-")
  eval(assignment[[1L]], env)
  expect_true(file.exists(env$streamline_figure_script))
  expect_identical(basename(env$streamline_figure_script), "fibrtraq-dtor985-plot-streamlines.R")
})

test_that("the installed dataset recipe builds its exploratory report with recorded methods", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("freesurferformats")
  script <- system.file("scripts", "fibrtraq-dtor985-start-to-finish.R", package = "FibrTraQ")
  expect_true(nzchar(script))
  assignment <- Filter(function(x) {
    is.call(x) && identical(x[[1L]], as.name("<-")) &&
      identical(x[[2L]], as.name("exploratory_report"))
  }, as.list(parse(script)))
  expect_length(assignment, 1L)

  fixture <- ftq_example_data(verify = TRUE)
  params <- fixture$parameters
  root <- tempfile("installed-recipe-exploratory-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  # Execute the recipe's exact report call on a real, tiny run. A nondefault
  # local-search radius proves that both public arguments reach the report.
  expect_warning(full_run <- ftq_run(list(
    dataset_root = fixture$root,
    subject = "recipe-check",
    description = "local-search",
    fnode = fixture$focal_mask,
    local = fixture$local_mask,
    tractogram = fixture$tractogram,
    ref_image = fixture$reference_image,
    targets = as.list(fixture$target_masks),
    coordinate = list(track_space = params$track_space,
      discretization = params$discretization,
      axis_flip = as.logical(params$axis_flip)),
    analysis = list(target_intersection_mode = "local_search",
      local_search_radius_mm = 1.5, endpoint_analysis = "none"),
    report = list(background_image = fixture$reference_image),
    outputs = list(output_dir = file.path(root, "run"),
      write_maps = FALSE, write_report = FALSE)
  ), backend = "r", verbose = FALSE, progress = FALSE),
  "does not declare NIfTI spatial units", fixed = TRUE)
  expect_true(full_run$validation$passed)

  env <- new.env(parent = environment())
  env$full_run <- full_run
  env$exploratory_output_dir <- file.path(root, "exploratory")
  env$report_prefix <- "recipe-check"
  env$OVERWRITE_OUTPUTS <- TRUE
  expect_no_warning(eval(assignment[[1L]], env))
  report <- env$exploratory_report
  expect_true(file.exists(report$html))
  expect_gt(file.info(report$html)$size, 0)
  provenance <- report$summary$provenance
  expect_null(report$summary$outputs$provenance_csv)
  expect_identical(provenance$target_intersection_mode, "local_search")
  expect_equal(provenance$local_search_radius_mm, 1.5)
  context <- attr(report$summary$provenance, "report_context")
  expect_identical(context$subject, "recipe-check")
  expect_identical(context$metadata_path, full_run$outputs$metadata_json)
  html <- paste(readLines(report$html, warn = FALSE), collapse = "\n")
  expect_match(html, "<b>Target assignment mode:</b> local_search", fixed = TRUE)
  expect_match(html, "<b>Local search radius (mm):</b> 1.5", fixed = TRUE)

  env$full_run$manifest$local_search_radius_mm <- 2
  expect_error(eval(assignment[[1L]], env), "conflicts with recorded result provenance")
})

test_that("new QC and repaired indices link current guides while retaining old guides", {
  root <- tempfile("guide revision # 100%20 ")
  viewer_dir <- file.path(root, "001", "focal", "reports")
  dir.create(viewer_dir, recursive = TRUE)
  old_guides <- file.path(rep(c(root, viewer_dir), each = 2L),
                         rep(c("OUTPUTS_README.html", "OUTPUTS_README.md"), 2L))
  for (path in old_guides) writeLines("Historical definitions retained.", path)
  old_hashes <- tools::md5sum(old_guides)
  fixture <- ftq_example_data()
  params <- fixture$parameters
  viewer <- ftq_interactive_qc_viewer(
    list(dataset_root = root, subject = "001", description = "guide-check",
         fnode = fixture$focal_mask, tractogram = fixture$tractogram,
         targets = as.list(fixture$target_masks),
         coordinate = list(track_space = params$track_space,
                           discretization = params$discretization,
                           axis_flip = as.logical(params$axis_flip)),
         outputs = list(output_dir = viewer_dir, write_maps = FALSE,
                        write_report = FALSE, overwrite = FALSE)),
    output_html = file.path(viewer_dir, "viewer.html"), sample_per_target = 1L
  )
  index <- data.frame(
    subject = "001", focal_region = "focal", run_status = "completed",
    validation_passed = TRUE, biological_qc_status = "unavailable",
    primary_report = viewer$html, exploratory_report = NA_character_,
    interactive_3d_qc = NA_character_, stringsAsFactors = FALSE
  )
  utils::write.csv(index, file.path(root, "fibrtraq-run-index.csv"), row.names = FALSE)
  repaired <- ftq_rebuild_report_index(root)
  for (html_path in c(viewer$html, repaired)) {
    text <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
    anchor <- regmatches(text, regexpr("href=['\"][^'\"]*OUTPUTS_README-[^'\"]+[.]html['\"]", text))
    expect_true(nzchar(anchor), info = html_path)
    href <- sub("^href=['\"](.*)['\"]$", "\\1", anchor)
    guide <- file.path(dirname(html_path), utils::URLdecode(href))
    expect_true(file.exists(guide))
    current <- paste(readLines(guide, warn = FALSE), collapse = "\n")
    expect_match(current, "FibrTraQ output guide", fixed = TRUE)
    expect_false(grepl("Historical definitions retained", current, fixed = TRUE))
  }
  expect_identical(tools::md5sum(old_guides), old_hashes)
})

test_that("primary navigation preserves content and travels with the complete bundle", {
  root <- tempfile("nav # 100%20 caf\u00e9 ")
  dir.create(file.path(root, "report"), recursive = TRUE)
  files <- c(main = file.path(root, "report", "main.html"),
             index = file.path(root, "index.html"),
             viewer = file.path(root, "report", "viewer # 100%20.html"))
  before <- "<!doctype html><h1>Subject 001 caf\u00e9</h1>"
  marker <- "<!-- FIBRTRAQ_NAVIGATION --><!-- /FIBRTRAQ_NAVIGATION -->"
  after <- "<p>Preserved report values &amp; labels.</p>\n"
  writeBin(charToRaw(enc2utf8(paste0(before, marker, after))), files[["main"]])
  writeLines("<!doctype html><title>Index</title>", files[["index"]])
  writeLines("<!doctype html><title>Viewer</title>", files[["viewer"]])
  links <- list("Results index" = files[["index"]],
                "3-D QC <viewer> & results" = files[["viewer"]],
                "Missing report" = file.path(root, "absent.html"))
  expect_true(ftq_report_add_navigation(files[["main"]], links))
  html <- rawToChar(readBin(files[["main"]], "raw", n = file.info(files[["main"]])$size))
  Encoding(html) <- "UTF-8"
  expect_match(html, "3-D QC &lt;viewer&gt; &amp; results", fixed = TRUE)
  expect_false(grepl("Missing report", html, fixed = TRUE))
  expect_identical(sub("<!-- FIBRTRAQ_NAVIGATION -->.*<!-- /FIBRTRAQ_NAVIGATION -->", marker, html),
                   paste0(before, marker, after))
  hash <- tools::md5sum(files[["main"]])
  expect_true(ftq_report_add_navigation(files[["main"]], links))
  expect_identical(tools::md5sum(files[["main"]]), hash)

  moved <- tempfile("moved nav # 100%20 ")
  dir.create(moved)
  expect_true(all(file.copy(list.files(root, full.names = TRUE), moved, recursive = TRUE)))
  new_main <- file.path(moved, "report", "main.html")
  hrefs <- regmatches(html, gregexpr("href='[^']+'", html))[[1L]]
  hrefs <- sub("^href='(.*)'$", "\\1", hrefs)
  expect_length(hrefs, 2L)
  expect_true(all(file.exists(file.path(dirname(new_main), utils::URLdecode(hrefs)))))
  expect_false(any(grepl("%5C|^file:|^/", hrefs)))

  historical <- file.path(root, "historical.html")
  writeLines("<!doctype html><p>No navigation marker</p>", historical)
  hash <- tools::md5sum(historical)
  expect_false(ftq_report_add_navigation(historical, links))
  expect_identical(tools::md5sum(historical), hash)
})
