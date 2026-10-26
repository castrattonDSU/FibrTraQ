# dTOR-985 FibrTraQ backend equivalence and speed benchmark.
#
# Run from RStudio with:
# source("inst/scripts/fibrtraq-dtor985-backend-benchmark.R")
#
# This script compares the reference R backend against the compiled C++ backend
# on the same dTOR-985 pilot and full Harvard-Oxford left-thalamus analyses.
# It records total elapsed runtime, internal stage timings, and exact output
# equivalence for the core biological result tables.

RUN_PILOT <- TRUE
RUN_FULL <- TRUE
BACKEND_ORDER <- c("r", "cpp")

# Keep maps and HTML off for benchmark fairness. This measures analysis time,
# not report rendering or NIfTI map-writing overhead.
WRITE_MAPS <- FALSE
WRITE_REPORT <- FALSE

OVERWRITE_OUTPUTS <- TRUE
SHOW_PROGRESS <- TRUE
PROGRESS_EVERY <- 100000L
STOP_ON_MISMATCH <- TRUE

repo_root <- Sys.getenv("FIBRTRAQ_REPO_ROOT", unset = normalizePath(getwd(), mustWork = FALSE))
data_root <- Sys.getenv("FIBRTRAQ_DATA_ROOT", unset = file.path(dirname(repo_root), "FibrTraQ-data"))
dtor_root <- Sys.getenv("FIBRTRAQ_DTOR_ROOT", unset = file.path(data_root, "public-validation", "dtor985"))

pilot_manifest <- file.path(dtor_root, "manifests", "dtor985-left-thalamus-harvardoxford-pilot.json")
full_manifest <- file.path(dtor_root, "manifests", "dtor985-left-thalamus-harvardoxford-full.json")

benchmark_root <- file.path(
  dtor_root,
  "outputs",
  "harvard-oxford-left-thalamus",
  "backend-benchmarks",
  format(Sys.time(), "%Y%m%d-%H%M%S")
)

announce <- function(text) {
  cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  cat(text, "\n")
  cat(paste(rep("=", 78), collapse = ""), "\n", sep = "")
}

require_existing_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path, call. = FALSE)
  invisible(path)
}

strip_ftq_attrs <- function(x) {
  attr(x, "fibrtraq_profile") <- NULL
  attr(x, "fibrtraq_backend") <- NULL
  x
}

profile_seconds <- function(profile, stage) {
  if (is.null(profile) || !nrow(profile)) return(NA_real_)
  value <- profile$seconds[profile$stage == stage]
  if (!length(value)) return(NA_real_)
  as.numeric(value[1])
}

target_columns <- function(result) {
  reserved <- c("fNode_Coords", "tracks", "local", "leave_no_target", "none", "cNode_Coords")
  setdiff(names(result), reserved)
}

summarize_run <- function(run) {
  result <- run$result
  targets <- target_columns(result)
  focal_voxel_streamline_associations <- sum(result$tracks)
  local_associations <- sum(result$local)
  leave_no_target_associations <- sum(result$leave_no_target)
  target_associations <- if (length(targets)) sum(unlist(result[targets], use.names = FALSE)) else 0
  data.frame(
    focal_voxels = nrow(result),
    focal_voxel_streamline_associations = focal_voxel_streamline_associations,
    local_associations = local_associations,
    leave_no_target_associations = leave_no_target_associations,
    leave_with_target_associations = focal_voxel_streamline_associations - local_associations - leave_no_target_associations,
    focal_voxel_streamline_target_associations = target_associations,
    validation_passed = isTRUE(run$validation$passed),
    stringsAsFactors = FALSE
  )
}

time_ftq_run <- function(label, manifest, backend) {
  output_dir <- file.path(benchmark_root, label, backend)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  announce(paste("Run", label, "with backend =", backend))
  started_at <- Sys.time()
  timing <- system.time({
    run <- ftq_run(
      manifest = manifest,
      output_dir = output_dir,
      overwrite = OVERWRITE_OUTPUTS,
      write_maps = WRITE_MAPS,
      write_report = WRITE_REPORT,
      verbose = TRUE,
      backend = backend,
      progress = SHOW_PROGRESS,
      progress_every = PROGRESS_EVERY,
      profile = TRUE
    )
  })
  completed_at <- Sys.time()

  summary <- summarize_run(run)
  profile <- run$profile
  if (!is.null(profile) && nrow(profile)) {
    profile$dataset <- label
    profile$requested_backend <- backend
    profile$output_dir <- output_dir
  }

  timing_row <- data.frame(
    dataset = label,
    requested_backend = backend,
    used_backend = as.character(run$backend),
    elapsed_seconds = as.numeric(timing[["elapsed"]]),
    user_seconds = as.numeric(timing[["user.self"]]),
    system_seconds = as.numeric(timing[["sys.self"]]),
    count_stage_seconds = profile_seconds(profile, "count_tracks"),
    prepare_stage_seconds = profile_seconds(profile, "prepare_rois_and_luts"),
    mapping_stage_seconds = profile_seconds(profile, "mapping_detection"),
    assemble_stage_seconds = profile_seconds(profile, "assemble_output"),
    started_at = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
    completed_at = format(completed_at, "%Y-%m-%dT%H:%M:%S%z"),
    result_csv = run$outputs$result_csv,
    rsc_csv = run$outputs$rsc_csv,
    fate_csv = run$outputs$fate_csv,
    validation_json = run$outputs$validation_json,
    output_dir = output_dir,
    stringsAsFactors = FALSE
  )
  list(
    run = run,
    timing = cbind(timing_row, summary),
    profile = profile
  )
}

compare_backend_pair <- function(label, runs) {
  r_run <- runs[["r"]]$run
  cpp_run <- runs[["cpp"]]$run

  result_identical <- identical(strip_ftq_attrs(r_run$result), strip_ftq_attrs(cpp_run$result))
  rsc_identical <- identical(r_run$rsc, cpp_run$rsc)

  r_fate <- utils::read.csv(r_run$outputs$fate_csv, stringsAsFactors = FALSE)
  cpp_fate <- utils::read.csv(cpp_run$outputs$fate_csv, stringsAsFactors = FALSE)
  fate_identical <- identical(r_fate, cpp_fate)

  comparison <- data.frame(
    dataset = label,
    result_identical = result_identical,
    rsc_identical = rsc_identical,
    fate_identical = fate_identical,
    all_core_outputs_identical = result_identical && rsc_identical && fate_identical,
    elapsed_speedup_cpp_vs_r = runs[["r"]]$timing$elapsed_seconds / runs[["cpp"]]$timing$elapsed_seconds,
    count_stage_speedup_cpp_vs_r = runs[["r"]]$timing$count_stage_seconds / runs[["cpp"]]$timing$count_stage_seconds,
    r_elapsed_seconds = runs[["r"]]$timing$elapsed_seconds,
    cpp_elapsed_seconds = runs[["cpp"]]$timing$elapsed_seconds,
    r_count_stage_seconds = runs[["r"]]$timing$count_stage_seconds,
    cpp_count_stage_seconds = runs[["cpp"]]$timing$count_stage_seconds,
    stringsAsFactors = FALSE
  )

  if (isTRUE(STOP_ON_MISMATCH) && !isTRUE(comparison$all_core_outputs_identical)) {
    stop("R and C++ outputs differed for ", label, ". See benchmark comparison table.", call. = FALSE)
  }

  comparison
}

run_benchmark_dataset <- function(label, manifest) {
  require_existing_file(manifest, paste(label, "manifest"))
  preflight <- ftq_preflight(manifest, overwrite = OVERWRITE_OUTPUTS, write_maps = WRITE_MAPS)
  if (!isTRUE(preflight$passed)) {
    print(preflight$checks)
    stop("Preflight failed for ", label, call. = FALSE)
  }

  runs <- list()
  timing_rows <- list()
  profile_rows <- list()

  for (backend in BACKEND_ORDER) {
    gc()
    backend_result <- time_ftq_run(label, manifest, backend)
    run <- backend_result$run
    row <- backend_result$timing
    runs[[backend]] <- list(run = run, timing = row)
    timing_rows[[backend]] <- row
    if (!is.null(backend_result$profile) && nrow(backend_result$profile)) {
      profile_rows[[backend]] <- transform(
        backend_result$profile,
        dataset = label,
        requested_backend = backend,
        output_dir = row$output_dir
      )
    }
  }

  comparison <- compare_backend_pair(label, runs)

  list(
    timing = do.call(rbind, timing_rows),
    profiles = if (length(profile_rows)) do.call(rbind, profile_rows) else data.frame(),
    comparison = comparison
  )
}

announce("Step 1: Load FibrTraQ")
if (requireNamespace("pkgload", quietly = TRUE) && file.exists(file.path(repo_root, "DESCRIPTION"))) {
  pkgload::load_all(repo_root, reset = TRUE, recompile = TRUE, export_all = FALSE, quiet = FALSE)
} else {
  library(FibrTraQ)
}

announce("Step 2: Verify benchmark inputs")
require_existing_file(pilot_manifest, "Pilot manifest")
require_existing_file(full_manifest, "Full manifest")
dir.create(benchmark_root, recursive = TRUE, showWarnings = FALSE)
cat("Benchmark output root:\n", benchmark_root, "\n", sep = "")

datasets <- list()
if (isTRUE(RUN_PILOT)) datasets$pilot <- pilot_manifest
if (isTRUE(RUN_FULL)) datasets$full <- full_manifest
if (!length(datasets)) stop("Set RUN_PILOT or RUN_FULL to TRUE.", call. = FALSE)

announce("Step 3: Run requested benchmarks")
results <- lapply(names(datasets), function(label) run_benchmark_dataset(label, datasets[[label]]))
names(results) <- names(datasets)

timing <- do.call(rbind, lapply(results, `[[`, "timing"))
profiles <- do.call(rbind, lapply(results, `[[`, "profiles"))
comparisons <- do.call(rbind, lapply(results, `[[`, "comparison"))

timing_csv <- file.path(benchmark_root, "dtor985_backend_timing.csv")
profile_csv <- file.path(benchmark_root, "dtor985_backend_profiles.csv")
comparison_csv <- file.path(benchmark_root, "dtor985_backend_equivalence.csv")

utils::write.csv(timing, timing_csv, row.names = FALSE)
utils::write.csv(profiles, profile_csv, row.names = FALSE)
utils::write.csv(comparisons, comparison_csv, row.names = FALSE)

announce("Step 4: Benchmark summary")
print(comparisons)
cat("\nTiming CSV:\n", timing_csv, "\n", sep = "")
cat("Profile CSV:\n", profile_csv, "\n", sep = "")
cat("Equivalence CSV:\n", comparison_csv, "\n", sep = "")

benchmark_outputs <- list(
  benchmark_root = benchmark_root,
  timing_csv = timing_csv,
  profile_csv = profile_csv,
  comparison_csv = comparison_csv,
  timing = timing,
  profiles = profiles,
  comparisons = comparisons
)

benchmark_outputs
