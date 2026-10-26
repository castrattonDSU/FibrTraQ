test_that("public validation registry distinguishes evidence from candidates", {
  registry_path <- system.file(
    "validation",
    "public-datasets.yml",
    package = "FibrTraQ"
  )
  expect_true(nzchar(registry_path))
  registry <- yaml::read_yaml(registry_path)
  datasets <- registry$datasets

  expect_true(length(datasets) > 0L)
  expect_equal(
    anyDuplicated(vapply(datasets, `[[`, character(1), "id")),
    0L
  )
  expect_false(any(vapply(
    datasets,
    function(dataset) "priority" %in% names(dataset),
    logical(1)
  )))

  allowed <- c(
    "completed_core_backend_benchmark",
    "template_only",
    "registry_only"
  )
  statuses <- vapply(
    datasets,
    `[[`,
    character(1),
    "evidence_status"
  )
  expect_true(all(statuses %in% allowed))
  expect_equal(
    statuses[vapply(datasets, `[[`, character(1), "id") == "dtor985"],
    "completed_core_backend_benchmark"
  )
})
