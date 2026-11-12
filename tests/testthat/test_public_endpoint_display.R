ftq_public_endpoint_plot_summary <- function() {
  data.frame(
    target = rep(c("connected", "empty"), 2L),
    tolerance = rep(c("exact", "dilated1"), each = 2L),
    target_intersection_streamlines = c(4, 0, 4, 0),
    focal_endpoint_streamlines = c(2, 2, 4, 2),
    paired_focal_target_endpoints = c(1, 0, 2, 0),
    shared_mask_ambiguous = c(1, 0, 1, 0),
    target_endpoint_only = c(0, 0, 0, 0),
    focal_endpoint_only = c(1, 0, 1, 0),
    neither_relevant_endpoint = c(1, 0, 0, 0),
    paired_percent_of_target_intersections = c(50, NA, 75, NA),
    paired_percent_of_focal_endpoint_streamlines = c(100, 0, 75, 0),
    stringsAsFactors = FALSE
  )
}

ftq_public_endpoint_plot_capture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$bars <- list()
  captured$labels <- character()
  captured$notes <- character()
  captured$target_labels <- list()
  captured$target_angles <- integer()
  original_axis <- graphics::axis
  original_text <- graphics::text
  original_mtext <- graphics::mtext
  exports <- list(
    `graphics::barplot` = function(height, ...) {
      captured$bars[[length(captured$bars) + 1L]] <- height
      args <- list(...)
      if (!is.null(args$names.arg)) {
        captured$target_labels[[length(captured$target_labels) + 1L]] <- args$names.arg
        captured$target_angles <- c(captured$target_angles, args$las)
      }
      draw_bars(height, ...)
    },
    `graphics::axis` = function(side, ...) {
      args <- list(...)
      if (side == 1 && !is.null(args$labels)) {
        captured$target_labels[[length(captured$target_labels) + 1L]] <- args$labels
        captured$target_angles <- c(captured$target_angles, args$las)
      }
      original_axis(side, ...)
    },
    `graphics::text` = function(x, y, labels, ...) {
      captured$labels <- c(captured$labels, rep_len(as.character(labels), max(length(x), length(y))))
      original_text(x, y, labels, ...)
    },
    `graphics::mtext` = function(text, ...) {
      captured$notes <- c(captured$notes, as.character(text))
      original_mtext(text, ...)
    }
  )
  # These plots supply numeric matrices. Instrument the real default renderer
  # as well, retaining coverage of its axis labels without changing graphics.
  draw_bars <- ftq_test_function_copy(graphics::barplot.default, bindings = list(
    axis = exports[["graphics::axis"]]
  ))
  bindings <- list(ftq_endpoint_plot_availability_note = ftq_test_function_copy(
    FibrTraQ:::ftq_endpoint_plot_availability_note, exports
  ))
  captured$classification <- ftq_test_function_copy(
    FibrTraQ:::ftq_plot_endpoint_classification, exports, bindings
  )
  captured$comparison <- ftq_test_function_copy(
    FibrTraQ:::ftq_plot_endpoint_comparison, exports, bindings
  )
  captured$penetration <- ftq_test_function_copy(
    FibrTraQ:::ftq_plot_endpoint_penetration, exports, bindings
  )
  captured
}

test_that("endpoint classification draws unavailable columns as NA, not zero percent", {
  captured <- ftq_public_endpoint_plot_capture()
  source <- ftq_public_endpoint_plot_summary()
  original <- source
  path <- tempfile("endpoint-classification-", fileext = ".png")
  on.exit(unlink(path), add = TRUE)
  expect_silent(captured$classification(source, path))
  expect_true(file.exists(path) && file.info(path)$size > 0)
  expect_equal(unname(captured$bars[[1L]][, 1L]), c(25, 25, 0, 25, 25))
  expect_true(all(is.na(captured$bars[[1L]][, 2L])))
  expect_true("NA" %in% captured$labels)
  expect_true(any(grepl("Unmarked zero-height bars mean 0%", captured$notes, fixed = TRUE)))
  expect_identical(source, original)
})

test_that("comparison retains geometric percentages and distinguishes each panel denominator", {
  captured <- ftq_public_endpoint_plot_capture()
  source <- ftq_public_endpoint_plot_summary()
  original <- source
  path <- tempfile("endpoint-comparison-", fileext = ".png")
  on.exit(unlink(path), add = TRUE)
  expect_silent(captured$comparison(source, path))
  expect_true(file.exists(path) && file.info(path)$size > 0)
  # Geometric percentages remain 50/75, rather than classification's 25/50.
  expect_equal(unname(captured$bars[[1L]][, 1L]), c(50, 75))
  expect_true(all(is.na(captured$bars[[1L]][, 2L])))
  expect_equal(unname(captured$bars[[2L]][, 1L]), c(100, 75))
  expect_equal(unname(captured$bars[[2L]][, 2L]), c(0, 0))
  expect_equal(sum(captured$labels == "NA"), 2L)
  expect_true(any(grepl("Geometric pairs include shared-mask membership", captured$notes, fixed = TRUE)))
  expect_identical(source, original)
})

test_that("all-undefined endpoint comparison remains a usable figure without warnings", {
  captured <- ftq_public_endpoint_plot_capture()
  source <- ftq_public_endpoint_plot_summary()
  source$target_intersection_streamlines <- 0
  source$focal_endpoint_streamlines <- 0
  # Legacy saved zeros must not imply that zero-denominator percentages exist.
  source$paired_percent_of_target_intersections <- 0
  source$paired_percent_of_focal_endpoint_streamlines <- 0
  original <- source
  path <- tempfile("endpoint-comparison-empty-", fileext = ".png")
  on.exit(unlink(path), add = TRUE)
  expect_silent(captured$comparison(source, path))
  expect_true(file.exists(path) && file.info(path)$size > 0)
  expect_true(all(vapply(captured$bars, function(x) all(is.na(x)), logical(1))))
  expect_equal(sum(captured$labels == "NA"), 8L)
  expect_identical(source, original)
})

test_that("penetration plots keep defined zero bins and mark absent intersection populations", {
  captured <- ftq_public_endpoint_plot_capture()
  source <- data.frame(
    target = c("connected", "empty"), target_intersection_streamlines = c(4, 0),
    one_voxel_touch_percent_of_target_intersections = c(25, 0),
    two_voxel_only_percent_of_target_intersections = c(0, 0),
    three_to_four_voxels_percent_of_target_intersections = c(25, 0),
    five_to_nine_voxels_percent_of_target_intersections = c(0, 0),
    ten_or_more_voxels_percent_of_target_intersections = c(50, 0)
  )
  original <- source
  path <- tempfile("endpoint-penetration-", fileext = ".png")
  on.exit(unlink(path), add = TRUE)
  expect_silent(captured$penetration(source, path))
  expect_true(file.exists(path) && file.info(path)$size > 0)
  expect_equal(unname(captured$bars[[1L]][, 1L]), c(25, 0, 25, 0, 50))
  expect_true(all(is.na(captured$bars[[1L]][, 2L])))
  expect_true("NA" %in% captured$labels)
  expect_true(any(grepl("not distance in millimetres", captured$notes, fixed = TRUE)))
  expect_identical(source, original)
})

test_that("all endpoint figures retain complete long target identities in readable labels", {
  captured <- ftq_public_endpoint_plot_capture()
  targets <- c("positive_x_cortical_GM", "negative_x_cortical_GM")
  summary <- ftq_public_endpoint_plot_summary()
  summary$target <- rep(targets, 2L)
  penetration <- data.frame(
    target = targets, target_intersection_streamlines = c(4, 4),
    one_voxel_touch_percent_of_target_intersections = c(25, 0),
    two_voxel_only_percent_of_target_intersections = c(0, 25),
    three_to_four_voxels_percent_of_target_intersections = c(25, 0),
    five_to_nine_voxels_percent_of_target_intersections = c(0, 25),
    ten_or_more_voxels_percent_of_target_intersections = c(50, 50)
  )
  original_summary <- summary
  original_penetration <- penetration
  paths <- vapply(seq_len(3L), function(i) tempfile(fileext = ".png"), character(1L))
  on.exit(unlink(paths), add = TRUE)
  expect_silent(captured$classification(summary, paths[[1L]]))
  expect_silent(captured$comparison(summary, paths[[2L]]))
  expect_silent(captured$penetration(penetration, paths[[3L]]))
  expect_true(length(captured$target_labels) >= 4L)
  for (labels in captured$target_labels) {
    expect_identical(unname(gsub("\n", "", labels, fixed = TRUE)), targets)
    expect_true(all(grepl("\n", labels, fixed = TRUE)))
  }
  expect_true(all(captured$target_angles == 1L))
  expect_identical(summary, original_summary)
  expect_identical(penetration, original_penetration)
})

test_that("an unbroken long target name cannot exceed the endpoint figure margins", {
  captured <- ftq_public_endpoint_plot_capture()
  source <- ftq_public_endpoint_plot_summary()
  targets <- c(paste0("positive", strrep("W", 192L)),
               paste0("negative", strrep("M", 192L)))
  source$target <- rep(targets, 2L)
  original <- source
  path <- tempfile("endpoint-long-token-", fileext = ".png")
  on.exit(unlink(path), add = TRUE)
  expect_silent(captured$comparison(source, path))
  expect_true(file.exists(path) && file.info(path)$size > 0)
  for (labels in captured$target_labels) {
    expect_identical(unname(gsub("\n", "", labels, fixed = TRUE)), targets)
    expect_true(all(nchar(unlist(strsplit(labels, "\n", fixed = TRUE))) <= 16L))
  }
  expect_identical(source, original)
})
