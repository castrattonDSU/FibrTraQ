test_that("sourcing the research archive registers its example print method", {
  # A caller outside the research environment must also find the method.
  method <- utils::getS3method("print", "ftq_example_result", envir = baseenv())
  expect_identical(method, research$print.ftq_example_result)
})
