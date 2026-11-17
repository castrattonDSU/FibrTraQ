test_that("the archived spacing audit measures anisotropic streamline segments", {
  tracks <- list(rbind(c(1, 1, 1), c(2, 1, 1), c(2, 3, 1)))
  audit <- ftq_spacing_chunk_cpp(tracks, diag(4), c(2, 3, 1),
                                bin_width_mm = 1, bin_width_vox = 1)
  expect_equal(c(audit$streamlines, audit$segments), c(1, 2))
  expect_equal(c(audit$min_vox, audit$max_vox), c(1, 2))
  expect_equal(c(audit$min_mm, audit$max_mm), c(2, 6))
  expect_equal(audit$hist_mm, c(0, 0, 1, 0, 0, 0, 1))
})
