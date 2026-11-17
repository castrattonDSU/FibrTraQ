#' Generate ISMRM 2015 validation manifests
#'
#' Read the ISMRM 2015 scoring-data segmentation configuration and create one
#' FibrTraQ manifest per bundle/sub-bundle with an available ground-truth `.trk`
#' file and configured head/tail endpoint masks.
#'
#' @param scoring_dir Path to the unpacked ISMRM 2015 scoring-data directory.
#'   This directory should contain `config_file_segmentation.json`, `bundles/`,
#'   `ROI/`, and usually `t1.nii.gz`.
#' @param output_dir Directory where manifest YAML files should be written.
#' @param bundles Optional character vector of bundle names to generate. If
#'   omitted, all configured bundles with matching `.trk` files are generated.
#' @param run_output_root Optional root directory for FibrTraQ outputs referenced
#'   inside the generated manifests. Defaults to
#'   `file.path(scoring_dir, "derivatives", derivative_name)`.
#' @param subject Subject label written into each manifest.
#' @param derivative_name Derivative/output name written into each manifest.
#' @param target_label Target ROI label used in the per-bundle manifest.
#' @param track_space Track coordinate-space setting written into each manifest.
#' @param discretization Discretization setting written into each manifest.
#' @param write_maps Whether generated manifests should request NIfTI maps.
#' @param write_report Whether generated manifests should request HTML reports.
#' @param overwrite Whether generated manifests should permit overwriting core
#'   output files by default.
#' @param include_background Whether to include `t1.nii.gz` as
#'   `report.background_image` when it exists.
#' @param force Whether existing manifest files may be overwritten.
#'
#' @return A data.frame describing generated manifests and input paths.
#' @export
ftq_make_ismrm2015_manifests <- function(scoring_dir,
                                         output_dir,
                                         bundles = NULL,
                                         run_output_root = NULL,
                                         subject = "sub-ismrm2015",
                                         derivative_name = "fibrtraq-ismrm2015",
                                         target_label = "tail",
                                         track_space = "auto",
                                         discretization = "any",
                                         write_maps = TRUE,
                                         write_report = TRUE,
                                         overwrite = FALSE,
                                         include_background = TRUE,
                                         force = FALSE) {
  scoring_dir <- normalizePath(path.expand(scoring_dir), mustWork = TRUE)
  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cfg_path <- file.path(scoring_dir, "config_file_segmentation.json")
  if (!file.exists(cfg_path)) {
    stop("ISMRM segmentation config does not exist: ", cfg_path)
  }
  cfg <- jsonlite::fromJSON(cfg_path, simplifyVector = FALSE)
  if (!is.list(cfg) || !length(cfg)) {
    stop("ISMRM segmentation config is empty or invalid: ", cfg_path)
  }

  all_bundles <- names(cfg)
  if (!is.null(bundles)) {
    bundles <- as.character(bundles)
    missing <- setdiff(bundles, all_bundles)
    if (length(missing)) {
      stop("Bundle(s) not present in ISMRM segmentation config: ",
           paste(missing, collapse = ", "))
    }
  } else {
    bundles <- all_bundles
  }

  if (is.null(run_output_root)) {
    run_output_root <- file.path(scoring_dir, "derivatives", derivative_name)
  }
  run_output_root <- normalizePath(path.expand(run_output_root), mustWork = FALSE)

  rows <- lapply(bundles, function(bundle) {
    item <- cfg[[bundle]]
    required <- c("head", "tail")
    missing_keys <- required[!required %in% names(item)]
    if (length(missing_keys)) {
      stop("Bundle ", bundle, " is missing config key(s): ",
           paste(missing_keys, collapse = ", "))
    }

    trk_rel <- ftq_ismrm2015_bundle_trk(scoring_dir, bundle)
    if (is.null(trk_rel)) {
      stop("No ground-truth .trk found for ISMRM bundle: ", bundle)
    }

    head_rel <- as.character(item$head)
    tail_rel <- as.character(item$tail)
    head_path <- file.path(scoring_dir, head_rel)
    tail_path <- file.path(scoring_dir, tail_rel)
    if (!file.exists(head_path)) stop("Missing head ROI for ", bundle, ": ", head_path)
    if (!file.exists(tail_path)) stop("Missing tail ROI for ", bundle, ": ", tail_path)

    bg_rel <- "t1.nii.gz"
    use_bg <- isTRUE(include_background) && file.exists(file.path(scoring_dir, bg_rel))
    desc <- paste0("ismrm2015-", bundle, "-groundtruth")
    manifest <- list(
      dataset_root = scoring_dir,
      subject = subject,
      description = desc,
      tractogram = trk_rel,
      fnode = head_rel,
      local = NULL,
      targets = stats::setNames(list(tail_rel), target_label),
      coordinate = list(
        track_space = track_space,
        discretization = discretization,
        axis_flip = "auto"
      ),
      outputs = list(
        derivative_name = derivative_name,
        output_dir = file.path(run_output_root, bundle),
        write_maps = isTRUE(write_maps),
        write_report = isTRUE(write_report),
        overwrite = isTRUE(overwrite)
      )
    )
    if (use_bg) {
      manifest$report <- list(background_image = bg_rel)
    }

    manifest_path <- file.path(
      output_dir,
      paste0("ismrm2015-", ftq_sanitize_label(bundle), ".yml")
    )
    if (file.exists(manifest_path) && !isTRUE(force)) {
      stop("Manifest already exists: ", manifest_path,
           ". Use force = TRUE to overwrite.")
    }
    yaml::write_yaml(manifest, manifest_path)

    data.frame(
      bundle = bundle,
      manifest_path = manifest_path,
      tractogram = file.path(scoring_dir, trk_rel),
      fnode = head_path,
      target = tail_path,
      output_dir = manifest$outputs$output_dir,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ftq_ismrm2015_bundle_trk <- function(scoring_dir, bundle) {
  candidates <- c(
    file.path("bundles", paste0(bundle, ".trk")),
    file.path("bundles", "sub_bundles", paste0(bundle, ".trk"))
  )
  found <- candidates[file.exists(file.path(scoring_dir, candidates))]
  if (!length(found)) return(NULL)
  found[1]
}

#' Generate ISMRM 2015 wrong-target control manifests
#'
#' Build negative-control manifests by pairing each ground-truth bundle's
#' tractogram and head endpoint with other bundles' tail endpoints. These runs
#' test whether FibrTraQ preserves endpoint specificity rather than counting any
#' focal streamline as a target connection.
#'
#' @inheritParams ftq_make_ismrm2015_manifests
#' @param target_bundles Optional character vector of bundles whose tail
#'   endpoints should be used as wrong targets. Defaults to all configured
#'   bundles except the source bundle for each manifest.
#' @param max_targets_per_bundle Optional integer cap used for quick pilot runs.
#'   If omitted, all eligible wrong targets are generated.
#' @param combine_targets If TRUE, generate one manifest per source bundle with
#'   all eligible wrong targets included as separate target columns. This gives
#'   the same target-specific counts with one tractogram pass per source bundle.
#'
#' @return A data.frame describing generated control manifests.
#' @export
ftq_make_ismrm2015_wrong_target_manifests <- function(scoring_dir,
                                                     output_dir,
                                                     bundles = NULL,
                                                     target_bundles = NULL,
                                                     max_targets_per_bundle = NULL,
                                                     run_output_root = NULL,
                                                     subject = "sub-ismrm2015",
                                                     derivative_name = "fibrtraq-ismrm2015-wrong-targets",
                                                     track_space = "auto",
                                                     discretization = "any",
                                                     write_maps = FALSE,
                                                     write_report = FALSE,
                                                     overwrite = FALSE,
                                                     include_background = FALSE,
                                                     combine_targets = FALSE,
                                                     force = FALSE) {
  ctx <- ftq_ismrm2015_context(scoring_dir, output_dir, run_output_root, derivative_name)
  cfg <- ctx$config
  all_bundles <- names(cfg)
  bundles <- ftq_ismrm2015_validate_bundles(bundles, all_bundles)
  target_bundles <- ftq_ismrm2015_validate_bundles(target_bundles, all_bundles)

  rows <- list()
  for (bundle in bundles) {
    src <- ftq_ismrm2015_bundle_record(ctx$scoring_dir, cfg, bundle)
    eligible <- setdiff(target_bundles, bundle)
    if (!is.null(max_targets_per_bundle)) eligible <- utils::head(eligible, as.integer(max_targets_per_bundle))
    if (isTRUE(combine_targets)) {
      targets <- lapply(eligible, function(target_bundle) {
        ftq_ismrm2015_bundle_record(ctx$scoring_dir, cfg, target_bundle)
      })
      target_paths <- vapply(targets, function(x) x$tail_rel, character(1))
      target_labels <- paste0("wrong_", eligible)
      desc <- paste0("ismrm2015-", bundle, "-wrong-target-matrix")
      manifest <- ftq_ismrm2015_manifest(
        scoring_dir = ctx$scoring_dir,
        subject = subject,
        description = desc,
        tractogram = src$trk_rel,
        fnode = src$head_rel,
        target_path = target_paths,
        target_label = target_labels,
        output_dir = file.path(ctx$run_output_root, bundle, "all-wrong-targets"),
        derivative_name = derivative_name,
        track_space = track_space,
        discretization = discretization,
        write_maps = write_maps,
        write_report = write_report,
        overwrite = overwrite,
        include_background = include_background
      )
      manifest_path <- file.path(
        ctx$output_dir,
        paste0("ismrm2015-", ftq_sanitize_label(bundle),
               "_wrong-target-matrix.yml")
      )
      ftq_ismrm2015_write_manifest(manifest, manifest_path, force)
      rows[[length(rows) + 1L]] <- data.frame(
        bundle = bundle,
        control = "wrong_target_matrix",
        target_bundle = paste(eligible, collapse = ";"),
        target_count = length(eligible),
        manifest_path = manifest_path,
        tractogram = src$trk_path,
        fnode = src$head_path,
        target = paste(vapply(targets, function(x) x$tail_path, character(1)), collapse = ";"),
        output_dir = manifest$outputs$output_dir,
        stringsAsFactors = FALSE
      )
      next
    }
    for (target_bundle in eligible) {
      tgt <- ftq_ismrm2015_bundle_record(ctx$scoring_dir, cfg, target_bundle)
      desc <- paste0("ismrm2015-", bundle, "-wrong-target-", target_bundle)
      manifest <- ftq_ismrm2015_manifest(
        scoring_dir = ctx$scoring_dir,
        subject = subject,
        description = desc,
        tractogram = src$trk_rel,
        fnode = src$head_rel,
        target_path = tgt$tail_rel,
        target_label = paste0("wrong_", target_bundle),
        output_dir = file.path(ctx$run_output_root, bundle, target_bundle),
        derivative_name = derivative_name,
        track_space = track_space,
        discretization = discretization,
        write_maps = write_maps,
        write_report = write_report,
        overwrite = overwrite,
        include_background = include_background
      )
      manifest_path <- file.path(
        ctx$output_dir,
        paste0("ismrm2015-", ftq_sanitize_label(bundle),
               "_wrong-target_", ftq_sanitize_label(target_bundle), ".yml")
      )
      ftq_ismrm2015_write_manifest(manifest, manifest_path, force)
      rows[[length(rows) + 1L]] <- data.frame(
        bundle = bundle,
        control = "wrong_target",
        target_bundle = target_bundle,
        manifest_path = manifest_path,
        tractogram = src$trk_path,
        fnode = src$head_path,
        target = tgt$tail_path,
        output_dir = manifest$outputs$output_dir,
        stringsAsFactors = FALSE
      )
    }
  }
  ftq_ismrm2015_rows(rows)
}

#' Generate ISMRM 2015 shifted-target control manifests
#'
#' Create non-wrapping shifted copies of each configured target endpoint mask and
#' build one control manifest per bundle. These controls test spatial
#' specificity: moving the target ROI should strongly reduce target associations.
#'
#' @inheritParams ftq_make_ismrm2015_manifests
#' @param shifted_roi_dir Directory where shifted target masks should be written.
#' @param shift Integer length-3 vector `c(dx, dy, dz)` in voxel units.
#'
#' @return A data.frame describing generated shifted-target manifests.
#' @export
ftq_make_ismrm2015_shifted_target_manifests <- function(scoring_dir,
                                                       output_dir,
                                                       shifted_roi_dir,
                                                       shift = c(8L, 0L, 0L),
                                                       bundles = NULL,
                                                       run_output_root = NULL,
                                                       subject = "sub-ismrm2015",
                                                       derivative_name = "fibrtraq-ismrm2015-shifted-targets",
                                                       target_label = "shifted_tail",
                                                       track_space = "auto",
                                                       discretization = "any",
                                                       write_maps = FALSE,
                                                       write_report = FALSE,
                                                       overwrite = FALSE,
                                                       include_background = FALSE,
                                                       force = FALSE) {
  ctx <- ftq_ismrm2015_context(scoring_dir, output_dir, run_output_root, derivative_name)
  cfg <- ctx$config
  all_bundles <- names(cfg)
  bundles <- ftq_ismrm2015_validate_bundles(bundles, all_bundles)
  shift <- ftq_ismrm2015_shift(shift)
  shifted_roi_dir <- normalizePath(path.expand(shifted_roi_dir), mustWork = FALSE)
  dir.create(shifted_roi_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  for (bundle in bundles) {
    rec <- ftq_ismrm2015_bundle_record(ctx$scoring_dir, cfg, bundle)
    shifted_path <- file.path(
      shifted_roi_dir,
      paste0(ftq_sanitize_label(bundle), "_tail_shift_",
             "dx", shift[1], "_dy", shift[2], "_dz", shift[3], ".nii.gz")
    )
    ftq_write_shifted_nifti(rec$tail_path, shifted_path, shift, force)
    desc <- paste0("ismrm2015-", bundle, "-shifted-target")
    manifest <- ftq_ismrm2015_manifest(
      scoring_dir = ctx$scoring_dir,
      subject = subject,
      description = desc,
      tractogram = rec$trk_rel,
      fnode = rec$head_rel,
      target_path = shifted_path,
      target_label = target_label,
      output_dir = file.path(ctx$run_output_root, bundle),
      derivative_name = derivative_name,
      track_space = track_space,
      discretization = discretization,
      write_maps = write_maps,
      write_report = write_report,
      overwrite = overwrite,
      include_background = include_background
    )
    manifest_path <- file.path(
      ctx$output_dir,
      paste0("ismrm2015-", ftq_sanitize_label(bundle), "_shifted-target.yml")
    )
    ftq_ismrm2015_write_manifest(manifest, manifest_path, force)
    rows[[length(rows) + 1L]] <- data.frame(
      bundle = bundle,
      control = "shifted_target",
      shift = paste(shift, collapse = ","),
      manifest_path = manifest_path,
      tractogram = rec$trk_path,
      fnode = rec$head_path,
      target = shifted_path,
      output_dir = manifest$outputs$output_dir,
      stringsAsFactors = FALSE
    )
  }
  ftq_ismrm2015_rows(rows)
}

#' Generate ISMRM 2015 submission sweep manifests
#'
#' Build manifests that run one or more submitted tractograms against every
#' configured ISMRM endpoint pair. This estimates which known bundle endpoints
#' are recovered by each submitted tractogram.
#'
#' @inheritParams ftq_make_ismrm2015_manifests
#' @param tractograms Character vector of submitted `.trk` files.
#'
#' @return A data.frame describing generated submission-sweep manifests.
#' @export
ftq_make_ismrm2015_submission_sweep_manifests <- function(scoring_dir,
                                                         tractograms,
                                                         output_dir,
                                                         bundles = NULL,
                                                         run_output_root = NULL,
                                                         subject = "sub-ismrm2015",
                                                         derivative_name = "fibrtraq-ismrm2015-submission-sweep",
                                                         target_label = "tail",
                                                         track_space = "auto",
                                                         discretization = "any",
                                                         write_maps = FALSE,
                                                         write_report = FALSE,
                                                         overwrite = FALSE,
                                                         include_background = FALSE,
                                                         force = FALSE) {
  ctx <- ftq_ismrm2015_context(scoring_dir, output_dir, run_output_root, derivative_name)
  cfg <- ctx$config
  bundles <- ftq_ismrm2015_validate_bundles(bundles, names(cfg))
  tractograms <- normalizePath(path.expand(tractograms), mustWork = TRUE)

  rows <- list()
  for (trk_path in tractograms) {
    submission <- tools::file_path_sans_ext(basename(trk_path))
    for (bundle in bundles) {
      rec <- ftq_ismrm2015_bundle_record(ctx$scoring_dir, cfg, bundle)
      desc <- paste0("ismrm2015-submission-", submission, "-", bundle)
      manifest <- ftq_ismrm2015_manifest(
        scoring_dir = ctx$scoring_dir,
        subject = subject,
        description = desc,
        tractogram = trk_path,
        fnode = rec$head_rel,
        target_path = rec$tail_rel,
        target_label = target_label,
        output_dir = file.path(ctx$run_output_root, submission, bundle),
        derivative_name = derivative_name,
        track_space = track_space,
        discretization = discretization,
        write_maps = write_maps,
        write_report = write_report,
        overwrite = overwrite,
        include_background = include_background
      )
      manifest_path <- file.path(
        ctx$output_dir,
        paste0("ismrm2015-submission-", ftq_sanitize_label(submission),
               "_", ftq_sanitize_label(bundle), ".yml")
      )
      ftq_ismrm2015_write_manifest(manifest, manifest_path, force)
      rows[[length(rows) + 1L]] <- data.frame(
        submission = submission,
        bundle = bundle,
        manifest_path = manifest_path,
        tractogram = trk_path,
        fnode = rec$head_path,
        target = rec$tail_path,
        output_dir = manifest$outputs$output_dir,
        stringsAsFactors = FALSE
      )
    }
  }
  ftq_ismrm2015_rows(rows)
}

ftq_ismrm2015_context <- function(scoring_dir, output_dir, run_output_root, derivative_name) {
  scoring_dir <- normalizePath(path.expand(scoring_dir), mustWork = TRUE)
  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cfg_path <- file.path(scoring_dir, "config_file_segmentation.json")
  if (!file.exists(cfg_path)) {
    stop("ISMRM segmentation config does not exist: ", cfg_path)
  }
  cfg <- jsonlite::fromJSON(cfg_path, simplifyVector = FALSE)
  if (!is.list(cfg) || !length(cfg)) {
    stop("ISMRM segmentation config is empty or invalid: ", cfg_path)
  }
  if (is.null(run_output_root)) {
    run_output_root <- file.path(scoring_dir, "derivatives", derivative_name)
  }
  list(
    scoring_dir = scoring_dir,
    output_dir = output_dir,
    run_output_root = normalizePath(path.expand(run_output_root), mustWork = FALSE),
    config = cfg
  )
}

ftq_ismrm2015_validate_bundles <- function(bundles, all_bundles) {
  if (is.null(bundles)) return(all_bundles)
  bundles <- as.character(bundles)
  missing <- setdiff(bundles, all_bundles)
  if (length(missing)) {
    stop("Bundle(s) not present in ISMRM segmentation config: ",
         paste(missing, collapse = ", "))
  }
  bundles
}

ftq_ismrm2015_bundle_record <- function(scoring_dir, cfg, bundle) {
  item <- cfg[[bundle]]
  required <- c("head", "tail")
  missing_keys <- required[!required %in% names(item)]
  if (length(missing_keys)) {
    stop("Bundle ", bundle, " is missing config key(s): ",
         paste(missing_keys, collapse = ", "))
  }
  trk_rel <- ftq_ismrm2015_bundle_trk(scoring_dir, bundle)
  if (is.null(trk_rel)) stop("No ground-truth .trk found for ISMRM bundle: ", bundle)
  head_rel <- as.character(item$head)
  tail_rel <- as.character(item$tail)
  out <- list(
    trk_rel = trk_rel,
    head_rel = head_rel,
    tail_rel = tail_rel,
    trk_path = file.path(scoring_dir, trk_rel),
    head_path = file.path(scoring_dir, head_rel),
    tail_path = file.path(scoring_dir, tail_rel)
  )
  if (!file.exists(out$head_path)) stop("Missing head ROI for ", bundle, ": ", out$head_path)
  if (!file.exists(out$tail_path)) stop("Missing tail ROI for ", bundle, ": ", out$tail_path)
  out
}

ftq_ismrm2015_manifest <- function(scoring_dir,
                                  subject,
                                  description,
                                  tractogram,
                                  fnode,
                                  target_path,
                                  target_label,
                                  output_dir,
                                  derivative_name,
                                  track_space,
                                  discretization,
                                  write_maps,
                                  write_report,
                                  overwrite,
                                  include_background) {
  target_path <- as.character(target_path)
  target_label <- as.character(target_label)
  if (length(target_path) != length(target_label)) {
    stop("target_path and target_label must have the same length.")
  }
  manifest <- list(
    dataset_root = scoring_dir,
    subject = subject,
    description = description,
    tractogram = tractogram,
    fnode = fnode,
    local = NULL,
    targets = stats::setNames(as.list(target_path), target_label),
    coordinate = list(
      track_space = track_space,
      discretization = discretization,
      axis_flip = "auto"
    ),
    outputs = list(
      derivative_name = derivative_name,
      output_dir = output_dir,
      write_maps = isTRUE(write_maps),
      write_report = isTRUE(write_report),
      overwrite = isTRUE(overwrite)
    )
  )
  bg_rel <- "t1.nii.gz"
  if (isTRUE(include_background) && file.exists(file.path(scoring_dir, bg_rel))) {
    manifest$report <- list(background_image = bg_rel)
  }
  manifest
}

ftq_ismrm2015_write_manifest <- function(manifest, manifest_path, force) {
  if (file.exists(manifest_path) && !isTRUE(force)) {
    stop("Manifest already exists: ", manifest_path,
         ". Use force = TRUE to overwrite.")
  }
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(manifest, manifest_path)
}

ftq_ismrm2015_rows <- function(rows) {
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ftq_ismrm2015_shift <- function(shift) {
  shift <- as.integer(shift)
  if (length(shift) != 3L || anyNA(shift)) {
    stop("shift must be an integer vector of length 3: c(dx, dy, dz).")
  }
  shift
}

ftq_write_shifted_nifti <- function(input_path, output_path, shift, force) {
  if (file.exists(output_path) && !isTRUE(force)) {
    stop("Shifted ROI already exists: ", output_path,
         ". Use force = TRUE to overwrite.")
  }
  template <- RNifti::readNifti(input_path, internal = FALSE)
  arr <- as.array(template)
  shifted <- ftq_shift_array(arr, shift)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  RNifti::writeNifti(shifted, output_path, template = template, datatype = "uint8")
  output_path
}

ftq_shift_array <- function(arr, shift) {
  dims <- dim(arr)[1:3]
  out <- array(0, dim = dims)
  src <- vector("list", 3L)
  dst <- vector("list", 3L)
  for (i in seq_len(3L)) {
    d <- dims[i]
    s <- shift[i]
    if (abs(s) >= d) return(out)
    if (s >= 0L) {
      src[[i]] <- seq_len(d - s)
      dst[[i]] <- (1L + s):d
    } else {
      src[[i]] <- (1L - s):d
      dst[[i]] <- seq_len(d + s)
    }
  }
  out[dst[[1]], dst[[2]], dst[[3]]] <- arr[src[[1]], src[[2]], src[[3]]]
  out
}
