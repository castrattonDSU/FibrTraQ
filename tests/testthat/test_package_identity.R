library(testthat)
library(FibrTraQ)

test_that("package identity handles source builds without remote fields", {
  package_path <- file.path(tempdir(), "identity-source", "FibrTraQ")
  identity <- FibrTraQ:::ftq_package_identity(
    description = list(
      Package = "FibrTraQ",
      Version = "0.3.0.9000"
    ),
    package_path = package_path
  )

  expect_named(identity, c(
    "name",
    "version",
    "remote_type",
    "remote_url",
    "remote_ref",
    "remote_sha",
    "built",
    "loaded_path"
  ))
  expect_identical(identity$name, "FibrTraQ")
  expect_identical(identity$version, "0.3.0.9000")
  expect_true(all(is.na(unlist(identity[c(
    "remote_type",
    "remote_url",
    "remote_ref",
    "remote_sha",
    "built"
  )]))))
  expect_identical(
    identity$loaded_path,
    normalizePath(package_path, mustWork = FALSE)
  )
})

test_that("package identity records GitHub fields without URL credentials", {
  identity <- FibrTraQ:::ftq_package_identity(
    description = list(
      Package = "FibrTraQ",
      Version = "0.4.0",
      RemoteType = "github",
      RemoteUrl = paste0(
        "https://credential-example-secret@github.com/",
        "castrattonDSU/FibrTraQ?credential=query-secret&download=1"
      ),
      RemoteRef = "v0.4.0",
      RemoteSha = "0123456789abcdef",
      Built = "R 4.5.3; aarch64-apple-darwin20; 2026-08-04"
    ),
    package_path = file.path(tempdir(), "identity-github", "FibrTraQ")
  )

  expect_identical(identity$remote_type, "github")
  expect_identical(
    identity$remote_url,
    "https://github.com/castrattonDSU/FibrTraQ"
  )
  expect_identical(identity$remote_ref, "v0.4.0")
  expect_identical(identity$remote_sha, "0123456789abcdef")
  expect_match(identity$built, "R 4.5.3", fixed = TRUE)
  expect_false(any(grepl(
    "example-secret|query-secret",
    unlist(identity),
    fixed = FALSE
  )))

  scp_identity <- FibrTraQ:::ftq_package_identity(
    description = list(
      Package = "FibrTraQ",
      Version = "0.4.0",
      RemoteUrl = paste0(
        "credential-scp-secret@github.com:",
        "castrattonDSU/FibrTraQ"
      )
    ),
    package_path = file.path(tempdir(), "identity-scp", "FibrTraQ")
  )
  expect_identical(
    scp_identity$remote_url,
    "github.com:castrattonDSU/FibrTraQ"
  )
  expect_false(grepl("scp-secret", scp_identity$remote_url, fixed = TRUE))

  remotes_github_identity <- FibrTraQ:::ftq_package_identity(
    description = list(
      Package = "FibrTraQ",
      Version = "0.4.0",
      RemoteType = "github",
      RemoteHost = "api.github.com",
      RemoteUsername = "castrattonDSU",
      RemoteRepo = "FibrTraQ",
      RemoteRef = "main",
      RemoteSha = "0123456789abcdef"
    ),
    package_path = file.path(tempdir(), "identity-remotes-github", "FibrTraQ")
  )
  expect_identical(
    remotes_github_identity$remote_url,
    "https://github.com/castrattonDSU/FibrTraQ"
  )

  legacy_github_identity <- FibrTraQ:::ftq_package_identity(
    description = list(
      Package = "FibrTraQ",
      Version = "0.4.0",
      GithubUsername = "castrattonDSU",
      GithubRepo = "FibrTraQ",
      GithubRef = "main",
      GithubSHA1 = "fedcba9876543210"
    ),
    package_path = file.path(tempdir(), "identity-legacy-github", "FibrTraQ")
  )
  expect_identical(legacy_github_identity$remote_type, "github")
  expect_identical(
    legacy_github_identity$remote_url,
    "https://github.com/castrattonDSU/FibrTraQ"
  )
  expect_identical(legacy_github_identity$remote_ref, "main")
  expect_identical(legacy_github_identity$remote_sha, "fedcba9876543210")

  malicious_fields <- list(
    list(
      RemoteHost = "credential@api.github.com",
      RemoteUsername = "castrattonDSU",
      RemoteRepo = "FibrTraQ"
    ),
    list(
      RemoteHost = "api.github.com",
      RemoteUsername = "credential@castrattonDSU",
      RemoteRepo = "FibrTraQ"
    ),
    list(
      RemoteHost = "api.github.com",
      RemoteUsername = "castrattonDSU",
      RemoteRepo = "FibrTraQ?credential=secret"
    )
  )
  for (fields in malicious_fields) {
    malicious_identity <- FibrTraQ:::ftq_package_identity(
      description = c(
        list(
          Package = "FibrTraQ",
          Version = "0.4.0",
          RemoteType = "github"
        ),
        fields
      ),
      package_path = file.path(tempdir(), "identity-malicious", "FibrTraQ")
    )
    expect_true(is.na(malicious_identity$remote_url))
    expect_false(any(grepl(
      "credential|secret",
      unlist(malicious_identity),
      ignore.case = TRUE
    )))
  }
})

test_that("written run metadata includes the loaded package identity", {
  output_dir <- tempfile("fibrtraq-package-identity-")
  result <- data.frame(
    fNode_Coords = "1_1_1",
    tracks = 1L,
    local = 0L,
    leave_no_target = 0L,
    T = 1L,
    none = 0L,
    stringsAsFactors = FALSE
  )
  rsc <- rsc_summary(result, percent = TRUE, warn_legacy = FALSE)

  outputs <- ftq_write_outputs(
    result = result,
    rsc = rsc,
    template_image = "unused.nii.gz",
    output_dir = output_dir,
    prefix = "sub-identity_desc-provenance",
    target_labels = "T",
    write_maps = FALSE,
    overwrite = TRUE
  )
  metadata <- jsonlite::fromJSON(
    outputs$metadata_json,
    simplifyVector = FALSE
  )
  expected <- FibrTraQ:::ftq_package_identity()

  expect_named(metadata$package, names(expected))
  expect_identical(metadata$package$name, expected$name)
  expect_identical(metadata$package$version, expected$version)
  expect_identical(metadata$package$loaded_path, expected$loaded_path)
})

test_that("regional run rows carry additive package identity provenance", {
  identity <- FibrTraQ:::ftq_package_identity(
    description = list(
      Package = "FibrTraQ",
      Version = "0.4.0",
      RemoteType = "github",
      RemoteUrl = "https://github.com/castrattonDSU/FibrTraQ",
      RemoteRef = "v0.4.0",
      RemoteSha = "0123456789abcdef",
      Built = "R 4.5.3; test-platform; 2026-08-04"
    ),
    package_path = file.path(tempdir(), "regional-identity", "FibrTraQ")
  )
  run <- list(
    package_identity = identity,
    manifest = list(
      target_intersection_mode = "local_search",
      local_search_radius_mm = 2,
      analysis = list()
    ),
    roi_summary = data.frame(
      processed_streamlines = 3,
      unique_focal_intersecting_streamlines = 2,
      selected_track_space = "voxel",
      selected_axis_flip = "FALSE,FALSE,FALSE",
      stringsAsFactors = FALSE
    ),
    validation = list(
      passed = TRUE,
      biological_qc_status = "passed"
    ),
    outputs = list(),
    report = NULL
  )

  row <- FibrTraQ:::ftq_regional_study_run_row(
    subject = "subject1",
    focal_region = "A",
    run = run,
    exploratory = NULL,
    viewer = NULL,
    started_at = "2026-08-04T12:00:00-0400",
    error = NULL
  )

  expect_identical(row$package_version, identity$version)
  expect_identical(row$package_remote_type, identity$remote_type)
  expect_identical(row$package_remote_url, identity$remote_url)
  expect_identical(row$package_remote_ref, identity$remote_ref)
  expect_identical(row$package_remote_sha, identity$remote_sha)
  expect_identical(row$package_built, identity$built)
  expect_identical(row$package_loaded_path, identity$loaded_path)
  expect_identical(row$target_intersection_mode, "local_search")
  expect_equal(row$local_search_radius_mm, 2)
})
