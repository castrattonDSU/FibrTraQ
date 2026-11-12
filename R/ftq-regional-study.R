#' Run one complete FibrTraQ analysis from files
#'
#' This is the simplest file-based FibrTraQ entry point. Supply one TrackVis or MRtrix
#' tractogram, one focal mask, one or more named target masks, and an output
#' directory. The function runs the standard analysis and writes the static
#' report, target report, endpoint tables, maps, interactive 3-D QC
#' viewer, and a report index.
#'
#' The function delegates to [ftq_run_regional_study()] with one subject and
#' one focal region. It does not implement separate counting semantics. Subject
#' identifiers remain character labels, including leading zeros. Subject and
#' region labels must produce distinct output folders after filename sanitization
#' and case folding, so the result bundle also works on Windows.
#'
#' Native `.tck` points are RAS+ world coordinates in millimeters. Their mapping
#' uses the NIfTI affine without an orientation search. The reference image, or
#' focal NIfTI when a reference is omitted, must declare spatial units and a
#' valid, invertible qform or sform; declared forms must agree when both are
#' present. Inputs are not automatically registered or resampled. Review
#' alignment independently before interpreting results. Counts are unweighted;
#' external SIFT2 weights are not applied.
#'
#' @param tractogram Path to one TrackVis `.trk` or MRtrix `.tck` file.
#'   TCK supports uncompressed, single-file Float32/Float64 data in either byte
#'   order. Compressed `.tck.gz` and external TCK payload files are unsupported.
#' @param focal_mask Path to one binary focal-region mask.
#' @param target_masks Named character vector or list of binary target-mask
#'   paths, for example `c(cerebrum = "...", pontine = "...")`.
#' @param output_dir Directory that will contain all derived outputs.
#' @param subject Subject label used in output paths and metadata. Supply the
#'   actual identifier explicitly, for example `subject = "S1c"` or `"001"`.
#'   If omitted, a message identifies the retained `"subject1"` default.
#' @param focal_label Label for the focal mask.
#' @param local_mask Optional binary local-region mask. When omitted, the focal
#'   mask is used as the package's local-region default.
#' @param reference_image Optional 3-D reference image, or a 4-D diffusion image
#'   accompanied by `bvals`. The reference is used for coordinate conversion
#'   when required and as the default report background. For TCK, the focal
#'   NIfTI supplies the coordinate reference when this is omitted. An independent
#'   anatomical reference remains strongly recommended for spatial review.
#' @param bvals Optional b-value file path or numeric vector. Required when
#'   `reference_image` or `background_image` is four-dimensional.
#' @param background_image Optional anatomical image for reports and the 3-D
#'   viewer. Defaults to `reference_image`.
#' @param endpoint_contexts Optional named character vector or list of trusted
#'   binary anatomical context-mask paths.
#' @param description Short description included in output prefixes.
#' @param low_b_max Maximum absolute b-value included when a 4-D image is
#'   reduced to a mean low-b reference.
#' @param reheader_image_to_mask_grid Whether a matching-dimension image with a
#'   different affine may be copied to the focal-mask grid without resampling.
#'   Leave this `FALSE` unless an independent geometry audit has established
#'   voxel-index correspondence and the source image header is non-authoritative.
#' @param track_space Coordinate-space setting passed to [ftq_run()]. Native
#'   TCK accepts `"auto"` or `"world"`; both use world coordinates deterministically.
#' @param discretization Discretization setting passed to [ftq_run()].
#' @param axis_flip Optional three-value logical axis flip. For TRK, when
#'   omitted, automatic mapping stops if its preferred candidate would reflect
#'   an axis. Verify alignment before explicitly setting `track_space` and
#'   `axis_flip`. Native TCK permits only `NULL` or three `FALSE` values.
#' @param focal_assignment_mode Focal assignment mode passed to [ftq_run()].
#' @param streamline_voxelization_mode Streamline voxelization mode passed to
#'   [ftq_run()].
#' @param endpoint_analysis Endpoint analysis mode passed to [ftq_run()].
#' @param track_reader Tractogram reader used by the run. `"auto"` streams both
#'   TRK and TCK by default; `"stream"` and `"memory"` may be selected explicitly.
#' @param track_chunk_size Streamline chunk size used by the streaming reader.
#' @param backend Counting backend passed to [ftq_run()].
#' @param write_maps Whether to write the standard NIfTI maps.
#' @param write_report Whether to write the standard static report.
#' @param write_exploratory_report Whether to write the target report.
#' @param write_qc_viewer Whether to write the interactive 3-D QC viewer.
#' @param sample_per_target Number of target-positive streamlines sampled per
#'   target in the interactive viewer.
#' @param overwrite Whether existing outputs may be replaced.
#' @param verbose Emit preparation and run messages.
#' @param progress Emit streamline-counting progress messages.
#' @param open Open the report index after successful completion.
#' @param assignment_mode Target assignment mode passed to [ftq_run()].
#'   Legacy and paper terminology aliases are accepted and normalized.
#' @param search_radius_mm Finite, strictly positive scanner-space search radius
#'   in millimeters used when `assignment_mode = "local_search"`.
#' @param report_slices Optional integer z-slice indices, using one-based
#'   focal-image array indexing, for report previews. `NULL` selects slices
#'   automatically. Every index must lie within each run's focal-image grid.
#'
#' @return An [ftq_run_regional_study()] result containing one run. The
#'   `report_index` member points to a local HTML page linking every report and
#'   the interactive viewer.
#' @export
ftq_run_files <- function(
    tractogram,
    focal_mask,
    target_masks,
    output_dir,
    subject = "subject1",
    focal_label = "focal",
    local_mask = NULL,
    reference_image = NULL,
    bvals = NULL,
    background_image = reference_image,
    endpoint_contexts = NULL,
    description = "analysis",
    low_b_max = 50,
    reheader_image_to_mask_grid = FALSE,
    track_space = "auto",
    discretization = "any",
    axis_flip = NULL,
    focal_assignment_mode = c("all_focal_voxels", "first_intersection"),
    streamline_voxelization_mode = c("sample_points", "segment_interpolation"),
    endpoint_analysis = c("exact_and_dilated1", "none"),
    track_reader = c("auto", "stream", "memory"),
    track_chunk_size = 50000L,
    backend = c("auto", "cpp", "r"),
    write_maps = TRUE,
    write_report = TRUE,
    write_exploratory_report = TRUE,
    write_qc_viewer = TRUE,
    sample_per_target = 25L,
    overwrite = FALSE,
    verbose = TRUE,
    progress = verbose,
    open = FALSE,
    assignment_mode = c(
      "full_streamline",
      "endpoint_voxels",
      "local_search"
    ),
    search_radius_mm = 2,
    report_slices = NULL) {
  target_intersection_mode <- ftq_normalize_target_intersection_mode(
    assignment_mode
  )
  local_search_radius_mm <- ftq_validate_local_search_radius(
    search_radius_mm
  )
  one_label <- function(x, argument) {
    x <- as.character(x)
    if (length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
      stop(argument, " must be one non-empty value.", call. = FALSE)
    }
    x
  }
  one_path <- function(x, argument) {
    x <- one_label(x, argument)
    normalizePath(path.expand(x), mustWork = TRUE)
  }

  if (missing(subject)) {
    message(
      "No subject identifier supplied: outputs will use 'subject1'. ",
      "Set subject explicitly, for example subject = 'S1c'."
    )
  }
  subject <- one_label(subject, "subject")
  focal_label <- one_label(focal_label, "focal_label")
  tractogram <- one_path(tractogram, "tractogram")
  focal_mask <- one_path(focal_mask, "focal_mask")
  target_masks <- ftq_regional_study_named_paths(
    target_masks,
    "target_masks",
    required_names = NULL
  )
  if (focal_label %in% names(target_masks)) {
    stop(
      "focal_label must differ from every target_masks label.",
      call. = FALSE
    )
  }

  subject_value <- function(value) {
    if (is.null(value)) return(NULL)
    stats::setNames(as.character(one_path(value, "image path")), subject)
  }
  bval_value <- if (is.null(bvals)) {
    NULL
  } else if (is.numeric(bvals)) {
    stats::setNames(list(as.numeric(bvals)), subject)
  } else {
    stats::setNames(one_path(bvals, "bvals"), subject)
  }
  local_value <- if (is.null(local_mask)) {
    NULL
  } else {
    stats::setNames(one_path(local_mask, "local_mask"), subject)
  }
  context_value <- if (is.null(endpoint_contexts)) {
    NULL
  } else {
    stats::setNames(list(endpoint_contexts), subject)
  }

  ftq_run_regional_study(
    tractograms = stats::setNames(tractogram, subject),
    region_masks = stats::setNames(
      list(c(stats::setNames(focal_mask, focal_label), target_masks)),
      subject
    ),
    output_dir = output_dir,
    reference_images = subject_value(reference_image),
    bvals = bval_value,
    background_images = subject_value(background_image),
    focal_regions = focal_label,
    local_masks = local_value,
    endpoint_contexts = context_value,
    description = description,
    low_b_max = low_b_max,
    reheader_images_to_region_grid = reheader_image_to_mask_grid,
    track_space = track_space,
    discretization = discretization,
    axis_flip = axis_flip,
    focal_assignment_mode = focal_assignment_mode,
    streamline_voxelization_mode = streamline_voxelization_mode,
    endpoint_analysis = endpoint_analysis,
    track_reader = track_reader,
    track_chunk_size = track_chunk_size,
    backend = backend,
    write_maps = write_maps,
    write_report = write_report,
    write_exploratory_report = write_exploratory_report,
    write_qc_viewer = write_qc_viewer,
    sample_per_target = sample_per_target,
    overwrite = overwrite,
    continue_on_error = FALSE,
    verbose = verbose,
    progress = progress,
    open = open,
    assignment_mode = target_intersection_mode,
    search_radius_mm = local_search_radius_mm,
    report_slices = report_slices
  )
}

#' Run a complete regional FibrTraQ study
#'
#' Run the standard FibrTraQ analysis for every requested focal region in every
#' subject, using all remaining regional masks as targets. The function keeps
#' study orchestration inside the package: it can prepare a 3-D mean low-b
#' reference from a 4-D diffusion image, build a same-grid regional union mask,
#' generate the primary and target reports, build the interactive 3-D QC
#' viewer, and write a study-level report index. Subject identifiers remain
#' character labels, including leading zeros. Subject and region labels must
#' produce distinct output folders after filename sanitization and case folding
#' so the result bundle also works on Windows.
#'
#' Native `.tck` inputs use RAS+ world coordinates in millimeters and require
#' valid NIfTI geometry. See [ftq_run_files()] for TCK format and spatial
#' requirements. The same reports, maps and viewer are generated for both formats.
#'
#' @param tractograms Named character vector or list of `.trk` or `.tck`
#'   tractogram paths, with one entry per subject.
#' @param region_masks Named list with one element per subject. Each subject
#'   element must be a named character vector or list of regional mask paths.
#'   Every subject must provide the same region labels.
#' @param output_dir Directory that will contain all derived study outputs.
#' @param reference_images Optional named character vector or list of reference
#'   image paths. Three-dimensional images are used directly. Four-dimensional
#'   images are reduced to a mean low-b image using `bvals`. For TCK, the focal
#'   NIfTI supplies the coordinate reference when an image is omitted. Independent
#'   anatomical references remain strongly recommended for spatial review.
#' @param bvals Optional named character vector/list of b-value file paths, or a
#'   named list of numeric b-value vectors. Required for each four-dimensional
#'   reference or background image.
#' @param background_images Optional named character vector/list of anatomical
#'   images for report and viewer context. When omitted, the prepared reference
#'   image is used.
#' @param focal_regions Region labels to analyze as focal regions. Defaults to
#'   every label in `region_masks`.
#' @param local_masks Either `"union"` to build a union of all supplied regional
#'   masks for each subject, `NULL` to use the focal mask as the package local
#'   default, or a named vector/list of subject-specific local-mask paths.
#' @param endpoint_contexts Optional named outer list with one element per
#'   subject. Each subject element is a named vector/list of binary anatomical
#'   context-mask paths.
#' @param description Short description included in output prefixes.
#' @param low_b_max Maximum absolute b-value included in a mean low-b reference.
#' @param reheader_images_to_region_grid Whether image voxel arrays with matching
#'   dimensions but different affines may be copied to the regional-mask grid
#'   without resampling. The default is `FALSE`. Set this only when a separate
#'   geometry audit has established voxel-index correspondence and the source
#'   image affine is known not to define the analysis grid. A named
#'   subject-specific logical vector/list is accepted.
#' @param track_space Coordinate-space setting passed to [ftq_run()]. Supply one
#'   value for all subjects or a named subject-specific vector/list. Native TCK
#'   accepts only `"auto"` or `"world"`; both select world coordinates directly.
#' @param discretization Discretization setting passed to [ftq_run()]. Supply
#'   one value/vector for all subjects or a named subject-specific list.
#' @param axis_flip Optional three-value logical axis flip. For TRK, when
#'   omitted, automatic mapping stops if its preferred candidate would reflect
#'   an axis. Verify alignment before explicitly setting `track_space` and
#'   `axis_flip`. A named subject-specific list is also accepted. Native TCK
#'   permits only `NULL` or three `FALSE` values.
#' @param focal_assignment_mode Focal assignment mode passed to [ftq_run()].
#' @param streamline_voxelization_mode Streamline voxelization mode passed to
#'   [ftq_run()].
#' @param endpoint_analysis Endpoint analysis mode passed to [ftq_run()].
#' @param track_reader Tractogram reader used by each run. `"auto"` streams both
#'   TRK and TCK by default; `"stream"` and `"memory"` may be selected explicitly.
#' @param track_chunk_size Streamline chunk size used by the streaming reader.
#' @param backend Counting backend passed to [ftq_run()].
#' @param write_maps Whether to write the standard NIfTI maps. Target-report
#'   figures are computed in memory and do not write exploratory NIfTI maps.
#' @param write_report Whether to write the standard static report.
#' @param write_exploratory_report Whether to write the target report.
#' @param write_qc_viewer Whether to write the interactive 3-D QC viewer.
#' @param sample_per_target Number of target-positive streamlines sampled per
#'   target in each interactive viewer.
#' @param overwrite Whether existing outputs may be replaced.
#' @param continue_on_error Whether to record an individual counting,
#'   report, or viewer failure and continue with the remaining subject/region
#'   combinations. Subject-level input preparation failures still stop the
#'   study before those runs begin. Completed analysis files and their separate
#'   validation/report statuses are saved to the index before a run error stops
#'   execution. A presentation failure does not discard scientific results.
#' @param verbose Emit study and run progress messages.
#' @param progress Emit streamline-counting progress messages.
#' @param open Open the study report index after successful completion.
#' @param assignment_mode Target assignment mode passed to [ftq_run()].
#'   Legacy and paper terminology aliases are accepted and normalized.
#' @param search_radius_mm Finite, strictly positive scanner-space search radius
#'   in millimeters used when `assignment_mode = "local_search"`.
#' @param report_slices Optional integer z-slice indices, using one-based
#'   focal-image array indexing, for report previews. `NULL` selects slices
#'   automatically. Every index must lie within each run's focal-image grid.
#'
#' @return An object containing the run index, prepared-input index, study
#'   report path, CSV index path, and output directory. The run index retains
#'   `run_status` and adds `analysis_status`, `validation_status`, individual
#'   presentation statuses and error messages. Biological QC is reported
#'   separately from computational validation.
#' @export
ftq_run_regional_study <- function(
    tractograms,
    region_masks,
    output_dir,
    reference_images = NULL,
    bvals = NULL,
    background_images = NULL,
    focal_regions = NULL,
    local_masks = "union",
    endpoint_contexts = NULL,
    description = "regional-study",
    low_b_max = 50,
    reheader_images_to_region_grid = FALSE,
    track_space = "auto",
    discretization = "any",
    axis_flip = NULL,
    focal_assignment_mode = c("all_focal_voxels", "first_intersection"),
    streamline_voxelization_mode = c("sample_points", "segment_interpolation"),
    endpoint_analysis = c("exact_and_dilated1", "none"),
    track_reader = c("auto", "stream", "memory"),
    track_chunk_size = 50000L,
    backend = c("auto", "cpp", "r"),
    write_maps = TRUE,
    write_report = TRUE,
    write_exploratory_report = TRUE,
    write_qc_viewer = TRUE,
    sample_per_target = 25L,
    overwrite = FALSE,
    continue_on_error = FALSE,
    verbose = TRUE,
    progress = verbose,
    open = FALSE,
    assignment_mode = c(
      "full_streamline",
      "endpoint_voxels",
      "local_search"
    ),
    search_radius_mm = 2,
    report_slices = NULL) {
  target_intersection_mode <- ftq_normalize_target_intersection_mode(
    assignment_mode
  )
  local_search_radius_mm <- ftq_validate_local_search_radius(
    search_radius_mm
  )
  focal_assignment_mode <- match.arg(focal_assignment_mode)
  streamline_voxelization_mode <- match.arg(streamline_voxelization_mode)
  endpoint_analysis <- match.arg(endpoint_analysis)
  track_reader <- match.arg(track_reader)
  backend <- match.arg(backend)
  report_slices <- ftq_report_slices(report_slices)

  low_b_max <- as.numeric(low_b_max)[1]
  track_chunk_size <- as.integer(track_chunk_size)[1]
  sample_per_target <- as.integer(sample_per_target)
  if (!is.finite(low_b_max) || low_b_max < 0) {
    stop("low_b_max must be one finite, non-negative number.", call. = FALSE)
  }
  if (is.na(track_chunk_size) || track_chunk_size < 1L) {
    stop("track_chunk_size must be a positive integer.", call. = FALSE)
  }
  if (!length(sample_per_target) || anyNA(sample_per_target) ||
      any(sample_per_target < 0L)) {
    stop("sample_per_target must contain non-negative integer values.", call. = FALSE)
  }
  if (isTRUE(write_qc_viewer) && !any(sample_per_target > 0L)) {
    stop(
      "sample_per_target must include a positive value when write_qc_viewer = TRUE.",
      call. = FALSE
    )
  }

  inputs <- ftq_regional_study_inputs(
    tractograms = tractograms,
    region_masks = region_masks,
    reference_images = reference_images,
    background_images = background_images,
    bvals = bvals,
    local_masks = local_masks,
    endpoint_contexts = endpoint_contexts
  )
  subjects <- inputs$subjects
  region_labels <- inputs$region_labels
  ftq_regional_study_check_labels(subjects, "subject")

  if (is.null(focal_regions)) {
    focal_regions <- region_labels
  } else {
    focal_regions <- unique(as.character(focal_regions))
    if (!length(focal_regions) || anyNA(focal_regions) ||
        any(!nzchar(focal_regions))) {
      stop("focal_regions must contain one or more non-empty region labels.",
           call. = FALSE)
    }
    missing_regions <- setdiff(focal_regions, region_labels)
    if (length(missing_regions)) {
      stop(
        "focal_regions are not present in region_masks: ",
        paste(missing_regions, collapse = ", "),
        call. = FALSE
      )
    }
  }

  ftq_regional_study_check_labels(focal_regions, "focal-region", reserved = "prepared")
  description <- ftq_sanitize_label(as.character(description)[1])
  if (is.na(description) || !nzchar(description)) {
    stop("description must be a non-empty value.", call. = FALSE)
  }
  if (identical(endpoint_analysis, "none") &&
      any(vapply(inputs$endpoint_contexts, length, integer(1)) > 0L)) {
    stop("endpoint_contexts require endpoint_analysis = 'exact_and_dilated1'.",
         call. = FALSE)
  }

  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  index_csv <- file.path(output_dir, "fibrtraq-run-index.csv")
  prepared_csv <- file.path(output_dir, "fibrtraq-prepared-inputs.csv")
  report_index <- file.path(output_dir, "fibrtraq-report-index.html")
  ftq_check_overwrite(c(index_csv, prepared_csv, report_index), overwrite)

  prepared_rows <- vector("list", length(subjects))
  names(prepared_rows) <- subjects
  run_rows <- vector("list", length(subjects) * length(focal_regions))
  run_number <- 0L
  total_runs <- length(run_rows)
  # Persist the available record after each stage, including before a default
  # stop. Existing scientific files are never reconstructed or replaced here.
  persist_index <- function() {
    prepared_index <- do.call(rbind, Filter(Negate(is.null), prepared_rows))
    if (is.null(prepared_index)) prepared_index <- data.frame(subject = character())
    run_index <- do.call(rbind, Filter(Negate(is.null), run_rows))
    if (is.null(run_index)) {
      run_index <- ftq_regional_study_run_row(
        subject = "", focal_region = "", run = NULL, exploratory = NULL,
        viewer = NULL, started_at = "", error = NULL
      )[FALSE, , drop = FALSE]
    }
    rownames(prepared_index) <- rownames(run_index) <- NULL
    utils::write.csv(prepared_index, prepared_csv, row.names = FALSE)
    utils::write.csv(run_index, index_csv, row.names = FALSE)
    ftq_write_regional_study_index(run_index, report_index, output_dir, description)
    invisible(list(prepared_inputs = prepared_index, index = run_index))
  }
  ftq_write_output_guide(output_dir)
  persist_index()

  for (subject in subjects) {
    if (isTRUE(verbose)) message("\nPreparing ", subject)
    masks <- inputs$region_masks[[subject]]
    template_path <- unname(masks[[1]])
    template <- RNifti::readNifti(template_path, internal = TRUE)
    ftq_report_slices(report_slices, dim(template)[3L])
    for (region in names(masks)) {
      mask_image <- RNifti::readNifti(masks[[region]], internal = TRUE)
      ftq_assert_same_grid(
        mask_image,
        template,
        context = paste0("region_masks[['", subject, "']][['", region, "']]"),
        template_context = paste0("region_masks[['", subject, "']][[1]]")
      )
    }

    subject_dir <- file.path(output_dir, ftq_sanitize_label(subject))
    prepared_dir <- file.path(subject_dir, "prepared")
    # Image and union-mask writers create this directory only when needed.

    reference_path <- NULL
    reference_source <- inputs$reference_images[[subject]]
    subject_reheader <- ftq_regional_study_subject_setting(
      reheader_images_to_region_grid, subject, subjects
    )
    subject_reheader <- as.logical(unlist(
      subject_reheader,
      use.names = FALSE
    ))
    if (length(subject_reheader) != 1L || is.na(subject_reheader)) {
      stop(
        "reheader_images_to_region_grid for subject '", subject,
        "' must be one logical value.",
        call. = FALSE
      )
    }
    if (!is.null(reference_source)) {
      reference_path <- ftq_regional_study_prepare_image(
        source_path = reference_source,
        bvals = inputs$bvals[[subject]],
        template_path = template_path,
        output_path = file.path(
          prepared_dir,
          paste0(ftq_sanitize_label(subject), "_prepared-reference.nii.gz")
        ),
        volume_csv = file.path(
          prepared_dir,
          paste0(ftq_sanitize_label(subject), "_prepared-reference-volumes.csv")
        ),
        low_b_max = low_b_max,
        reheader_to_template = subject_reheader,
        overwrite = overwrite,
        context = paste0("reference_images[['", subject, "']]")
      )
    }

    background_source <- inputs$background_images[[subject]]
    background_path <- if (is.null(background_source)) {
      reference_path
    } else if (!is.null(reference_source) &&
               identical(background_source, reference_source)) {
      reference_path
    } else {
      ftq_regional_study_prepare_image(
        source_path = background_source,
        bvals = inputs$bvals[[subject]],
        template_path = template_path,
        output_path = file.path(
          prepared_dir,
          paste0(ftq_sanitize_label(subject), "_prepared-background.nii.gz")
        ),
        volume_csv = file.path(
          prepared_dir,
          paste0(ftq_sanitize_label(subject), "_prepared-background-volumes.csv")
        ),
        low_b_max = low_b_max,
        reheader_to_template = subject_reheader,
        overwrite = overwrite,
        context = paste0("background_images[['", subject, "']]")
      )
    }

    local_path <- NULL
    if (identical(inputs$local_mode, "union")) {
      local_path <- file.path(
        prepared_dir,
        paste0(ftq_sanitize_label(subject), "_regional-union.nii.gz")
      )
      union_inputs_path <- file.path(
        prepared_dir,
        paste0(ftq_sanitize_label(subject), "_regional-union-inputs.csv")
      )
      union_overlaps_path <- file.path(
        prepared_dir,
        paste0(ftq_sanitize_label(subject), "_regional-union-overlaps.csv")
      )
      ftq_check_overwrite(
        c(local_path, union_inputs_path, union_overlaps_path),
        overwrite
      )
      union_result <- ftq_make_union_mask(
        masks = as.list(masks),
        template_image = template_path,
        output_path = local_path,
        label = "regional_union",
        overwrite = overwrite
      )
      utils::write.csv(
        union_result$inputs,
        union_inputs_path,
        row.names = FALSE
      )
      utils::write.csv(
        union_result$pairwise_overlaps,
        union_overlaps_path,
        row.names = FALSE
      )
    } else if (identical(inputs$local_mode, "supplied")) {
      local_path <- inputs$local_masks[[subject]]
      local_image <- RNifti::readNifti(local_path, internal = TRUE)
      ftq_assert_same_grid(
        local_image,
        template,
        context = paste0("local_masks[['", subject, "']]"),
        template_context = paste0("region_masks[['", subject, "']][[1]]")
      )
    }

    prepared_rows[[subject]] <- data.frame(
      subject = subject,
      tractogram = inputs$tractograms[[subject]],
      reference_source = ftq_regional_study_scalar(reference_source),
      reference_source_md5 = ftq_regional_study_checksum(reference_source),
      bvals_source = ftq_regional_study_bval_source(inputs$bvals[[subject]]),
      bvals_source_md5 = ftq_regional_study_checksum(
        ftq_regional_study_bval_source(inputs$bvals[[subject]])
      ),
      reheadered_to_region_grid = subject_reheader,
      reference_image = ftq_regional_study_scalar(reference_path),
      background_source = ftq_regional_study_scalar(background_source),
      background_source_md5 = ftq_regional_study_checksum(background_source),
      background_image = ftq_regional_study_scalar(background_path),
      local_mask = ftq_regional_study_scalar(local_path),
      stringsAsFactors = FALSE
    )

    persist_index()

    for (focal_region in focal_regions) {
      run_number <- run_number + 1L
      if (isTRUE(verbose)) {
        message(
          "\n[", run_number, "/", total_runs, "] ",
          subject, " / ", focal_region
        )
      }

      run_output_dir <- file.path(subject_dir, ftq_sanitize_label(focal_region))
      subject_track_space <- ftq_regional_study_subject_setting(
        track_space, subject, subjects
      )
      subject_discretization <- ftq_regional_study_subject_setting(
        discretization, subject, subjects
      )
      subject_axis_flip <- ftq_regional_study_subject_setting(
        axis_flip, subject, subjects
      )
      if (!is.null(subject_axis_flip)) {
        subject_axis_flip <- as.logical(unlist(subject_axis_flip, use.names = FALSE))
        if (length(subject_axis_flip) != 3L || anyNA(subject_axis_flip)) {
          stop(
            "axis_flip for subject '", subject,
            "' must contain exactly three logical values.",
            call. = FALSE
          )
        }
      }

      manifest <- list(
        dataset_root = output_dir,
        subject = subject,
        focal_region = focal_region,
        description = paste0(
          description, "-", ftq_sanitize_label(focal_region),
          "-regional-connectivity"
        ),
        tractogram = inputs$tractograms[[subject]],
        fnode = masks[[focal_region]],
        local = local_path,
        ref_image = reference_path,
        targets = as.list(masks[setdiff(region_labels, focal_region)]),
        endpoint_contexts = inputs$endpoint_contexts[[subject]],
        coordinate = list(
          track_space = subject_track_space,
          discretization = subject_discretization,
          axis_flip = if (is.null(subject_axis_flip)) "auto" else
            as.list(subject_axis_flip)
        ),
        analysis = list(
          focal_assignment_mode = focal_assignment_mode,
          streamline_voxelization_mode = streamline_voxelization_mode,
          target_intersection_mode = target_intersection_mode,
          local_search_radius_mm = local_search_radius_mm,
          endpoint_analysis = endpoint_analysis
        ),
        execution = list(
          track_reader = track_reader,
          track_chunk_size = track_chunk_size
        ),
        outputs = list(
          output_dir = run_output_dir,
          write_maps = isTRUE(write_maps),
          write_report = isTRUE(write_report),
          overwrite = isTRUE(overwrite)
        ),
        report = list(
          background_image = background_path,
          context_source = if (is.null(background_path)) "auto" else "anatomy"
        )
      )

      started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      run <- exploratory <- viewer <- NULL
      stage_errors <- list(analysis = NULL, validation = NULL, primary_report = NULL,
                           exploratory_report = NULL, interactive_3d_qc = NULL)
      capture_stage <- function(stage, expression) {
        tryCatch(force(expression), error = function(e) {
          stage_errors[[stage]] <<- conditionMessage(e)
          NULL
        })
      }
      run <- capture_stage("analysis", ftq_run(
        manifest = manifest,
        overwrite = overwrite,
        write_maps = write_maps,
        write_report = write_report,
        verbose = verbose,
        backend = backend,
        progress = progress,
        report_slices = report_slices,
        report_error_action = "return"
      ))
      if (!is.null(run)) {
        if (!isTRUE(run$validation$passed)) {
          stage_errors$validation <-
            "The analysis completed but computational validation did not pass."
        }
        if (!is.null(run$report_error) && nzchar(run$report_error)) {
          stage_errors$primary_report <- run$report_error
        }
        # Independent presentation stages can still provide useful access even
        # when another report fails. Their status never changes validation.
        if (isTRUE(write_exploratory_report)) {
          exploratory <- capture_stage("exploratory_report", ftq_exploratory_report(
            result_csv = run$outputs$result_csv,
            output_dir = run_output_dir,
            template_image = masks[[focal_region]],
            prefix = run$manifest$prefix,
            target_labels = run$manifest$target_labels,
            background_image = background_path,
            tractogram_context = run$outputs$context,
            overwrite = overwrite,
            assignment_mode = run$manifest$target_intersection_mode,
            search_radius_mm = run$manifest$local_search_radius_mm,
            report_slices = report_slices
          ))
        }
        if (isTRUE(write_qc_viewer)) {
          viewer <- capture_stage("interactive_3d_qc", ftq_interactive_qc_viewer(
            run,
            sample_per_target = sample_per_target,
            overwrite = overwrite,
            open = FALSE
          ))
        }
      }
      if (!is.null(run$report$html) && file.exists(run$report$html)) {
        capture_stage("primary_report", ftq_report_add_navigation(
          run$report$html,
          list("Results index" = report_index,
               "Target report" = exploratory$html,
               "3-D QC viewer" = viewer$html)
        ))
      }
      failed_stages <- names(stage_errors)[lengths(stage_errors) > 0L]
      stage_labels <- c(analysis = "Analysis", validation = "Computational validation",
                        primary_report = "Main report", exploratory_report = "Target report",
                        interactive_3d_qc = "3D QC viewer")
      run_error <- if (length(failed_stages)) paste(vapply(
        failed_stages, function(stage) paste0(stage_labels[[stage]], ": ", stage_errors[[stage]]),
        character(1)
      ), collapse = "\n") else NULL

      run_rows[[run_number]] <- ftq_regional_study_run_row(
        subject = subject,
        focal_region = focal_region,
        run = run,
        exploratory = exploratory,
        viewer = viewer,
        started_at = started_at,
        error = run_error,
        manifest = manifest,
        stage_errors = stage_errors,
        requested = c(primary_report = isTRUE(write_report),
                      exploratory_report = isTRUE(write_exploratory_report),
                      interactive_3d_qc = isTRUE(write_qc_viewer))
      )
      persist_index()
      if (!is.null(run_error) && !isTRUE(continue_on_error)) {
        stop(
          "Regional study failed for subject '", subject,
          "', focal region '", focal_region, "':\n", run_error,
          "\nCompleted files and stage statuses are available in: ", report_index,
          call. = FALSE
        )
      }

      rm(run, exploratory, viewer)
      invisible(gc(verbose = FALSE))
    }
  }

  saved <- persist_index()
  prepared_index <- saved$prepared_inputs
  run_index <- saved$index

  out <- structure(
    list(
      index = run_index,
      prepared_inputs = prepared_index,
      index_csv = index_csv,
      prepared_inputs_csv = prepared_csv,
      report_index = report_index,
      output_dir = output_dir
    ),
    class = "ftq_regional_study_result"
  )
  if (isTRUE(open)) ftq_qc_open_output(report_index)
  out
}

#' @export
print.ftq_regional_study_result <- function(x, ...) {
  completed <- sum(x$index$run_status == "completed", na.rm = TRUE)
  failed <- sum(x$index$run_status == "failed", na.rm = TRUE)
  analyses <- if ("analysis_status" %in% names(x$index)) {
    sum(x$index$analysis_status == "completed", na.rm = TRUE)
  } else completed
  cat(
    "FibrTraQ results\n",
    "  Analyses saved: ", analyses, "\n",
    "  Completed runs: ", completed, "\n",
    "  Failed runs:    ", failed, "\n",
    "  Output:         ", x$output_dir, "\n",
    "  Report index:   ", x$report_index, "\n",
    sep = ""
  )
  invisible(x)
}

ftq_regional_study_inputs <- function(tractograms,
                                      region_masks,
                                      reference_images,
                                      background_images,
                                      bvals,
                                      local_masks,
                                      endpoint_contexts) {
  tractograms <- ftq_regional_study_named_paths(
    tractograms, "tractograms", required_names = NULL
  )
  subjects <- names(tractograms)
  if (!length(subjects)) {
    stop("tractograms must contain at least one subject.", call. = FALSE)
  }

  if (!is.list(region_masks) || is.null(names(region_masks)) || anyNA(names(region_masks)) ||
      any(!nzchar(names(region_masks))) || anyDuplicated(names(region_masks))) {
    stop(
      "region_masks must be a named list with one element per subject.",
      call. = FALSE
    )
  }
  ftq_regional_study_require_subjects(names(region_masks), subjects, "region_masks")
  region_masks <- region_masks[subjects]
  region_masks <- lapply(subjects, function(subject) {
    ftq_regional_study_named_paths(
      region_masks[[subject]],
      paste0("region_masks[['", subject, "']]"),
      required_names = NULL
    )
  })
  names(region_masks) <- subjects
  if (any(vapply(region_masks, length, integer(1)) < 2L)) {
    stop("Each subject in region_masks must contain at least two regions.",
         call. = FALSE)
  }
  region_labels <- names(region_masks[[1]])
  inconsistent <- subjects[!vapply(
    region_masks,
    function(x) setequal(names(x), region_labels),
    logical(1)
  )]
  if (length(inconsistent)) {
    stop(
      "Every subject must provide the same region labels. Mismatch: ",
      paste(inconsistent, collapse = ", "),
      call. = FALSE
    )
  }
  region_masks <- lapply(region_masks, function(x) x[region_labels])

  reference_images <- ftq_regional_study_optional_subject_paths(
    reference_images, subjects, "reference_images"
  )
  background_images <- ftq_regional_study_optional_subject_paths(
    background_images, subjects, "background_images"
  )
  bvals <- ftq_regional_study_bvals(bvals, subjects)

  if (is.null(local_masks)) {
    local_mode <- "none"
    local_masks <- stats::setNames(vector("list", length(subjects)), subjects)
  } else if (is.character(local_masks) && length(local_masks) == 1L &&
             identical(tolower(local_masks), "union")) {
    local_mode <- "union"
    local_masks <- stats::setNames(vector("list", length(subjects)), subjects)
  } else {
    local_mode <- "supplied"
    local_paths <- ftq_regional_study_named_paths(
      local_masks, "local_masks", required_names = subjects
    )
    local_masks <- as.list(local_paths[subjects])
  }

  endpoint_contexts <- ftq_regional_study_contexts(endpoint_contexts, subjects)
  list(
    subjects = subjects,
    region_labels = region_labels,
    tractograms = tractograms,
    region_masks = region_masks,
    reference_images = reference_images,
    background_images = background_images,
    bvals = bvals,
    local_mode = local_mode,
    local_masks = local_masks,
    endpoint_contexts = endpoint_contexts
  )
}

ftq_regional_study_named_paths <- function(x, argument, required_names = NULL) {
  if (is.list(x)) x <- unlist(x, use.names = TRUE)
  if (!is.character(x) || !length(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x))) {
    stop(
      argument,
      " must be a non-empty named character vector or named list of paths.",
      call. = FALSE
    )
  }
  if (!is.null(required_names)) {
    ftq_regional_study_require_subjects(names(x), required_names, argument)
    x <- x[required_names]
  }
  paths <- vapply(x, function(path) {
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
      stop(argument, " contains an empty or invalid path.", call. = FALSE)
    }
    normalizePath(path.expand(path), mustWork = TRUE)
  }, character(1))
  stats::setNames(unname(paths), names(x))
}

ftq_regional_study_optional_subject_paths <- function(x, subjects, argument) {
  out <- stats::setNames(vector("list", length(subjects)), subjects)
  if (is.null(x)) return(out)
  paths <- ftq_regional_study_named_paths(
    x, argument, required_names = subjects
  )
  as.list(paths[subjects])
}

ftq_regional_study_require_subjects <- function(actual, expected, argument) {
  missing <- setdiff(expected, actual)
  extra <- setdiff(actual, expected)
  if (length(missing) || length(extra)) {
    detail <- c(
      if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(extra)) paste0("unexpected: ", paste(extra, collapse = ", "))
    )
    stop(
      argument, " subject names must match tractograms (",
      paste(detail, collapse = "; "), ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

ftq_regional_study_bvals <- function(x, subjects) {
  out <- stats::setNames(vector("list", length(subjects)), subjects)
  if (is.null(x)) return(out)
  if (is.character(x) && !is.list(x)) {
    paths <- ftq_regional_study_named_paths(
      x, "bvals", required_names = subjects
    )
    return(as.list(paths[subjects]))
  }
  if (!is.list(x) || is.null(names(x)) || any(!nzchar(names(x))) ||
      anyDuplicated(names(x))) {
    stop(
      "bvals must be a named vector of file paths or a named list of paths/numeric vectors.",
      call. = FALSE
    )
  }
  ftq_regional_study_require_subjects(names(x), subjects, "bvals")
  for (subject in subjects) {
    value <- x[[subject]]
    if (is.character(value) && length(value) == 1L) {
      out[[subject]] <- normalizePath(path.expand(value), mustWork = TRUE)
    } else if (is.numeric(value)) {
      if (!length(value) || any(!is.finite(value))) {
        stop("bvals[['", subject, "']] contains invalid values.", call. = FALSE)
      }
      out[[subject]] <- as.numeric(value)
    } else {
      stop(
        "bvals[['", subject,
        "']] must be a file path or numeric vector.",
        call. = FALSE
      )
    }
  }
  out
}

ftq_regional_study_contexts <- function(x, subjects) {
  out <- stats::setNames(vector("list", length(subjects)), subjects)
  if (is.null(x)) return(out)
  if (!is.list(x) || is.null(names(x)) || any(!nzchar(names(x))) ||
      anyDuplicated(names(x))) {
    stop("endpoint_contexts must be a named list keyed by subject.",
         call. = FALSE)
  }
  unknown <- setdiff(names(x), subjects)
  if (length(unknown)) {
    stop(
      "endpoint_contexts contains unknown subject(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  for (subject in intersect(subjects, names(x))) {
    if (is.null(x[[subject]]) || !length(x[[subject]])) next
    out[[subject]] <- as.list(ftq_regional_study_named_paths(
      x[[subject]],
      paste0("endpoint_contexts[['", subject, "']]"),
      required_names = NULL
    ))
  }
  out
}

ftq_regional_study_prepare_image <- function(source_path,
                                             bvals,
                                             template_path,
                                             output_path,
                                             volume_csv,
                                             low_b_max,
                                             reheader_to_template,
                                             overwrite,
                                             context) {
  source <- RNifti::readNifti(source_path, internal = TRUE)
  template <- RNifti::readNifti(template_path, internal = TRUE)
  source_units <- ftq_local_search_unit_info(source)
  template_units <- ftq_local_search_unit_info(template)
  write_template <- template
  if (!template_units$spatial_units_known && source_units$spatial_units_known) {
    write_template <- RNifti::readNifti(template_path, internal = FALSE)
    RNifti::pixunits(write_template) <- source_units$spatial_units
  }
  dimensions_match <- identical(
    as.integer(dim(source)[1:3]),
    as.integer(dim(template)[1:3])
  )
  grid_matches <- dimensions_match && ftq_affine_equal(source, template)
  if (!dimensions_match) {
    ftq_assert_same_grid(
      source,
      template,
      context = context,
      template_context = "regional-mask template"
    )
  }
  if (!grid_matches && !isTRUE(reheader_to_template)) {
    stop(
      "dimensions match but affine/grid differs between ", context,
      " and the regional-mask template. Correct or regrid the input using ",
      "trusted geometry. If an independent audit proves that the voxel arrays ",
      "already correspond and only the source header is non-authoritative, ",
      "set reheader_images_to_region_grid = TRUE to create an explicit ",
      "non-resampled derivative.",
      call. = FALSE
    )
  }

  source_dims <- dim(source)
  nonspatial_dims <- if (length(source_dims) > 3L) {
    as.integer(source_dims[-seq_len(3L)])
  } else {
    integer()
  }
  volume_count <- if (length(nonspatial_dims)) prod(nonspatial_dims) else 1L
  if (volume_count == 1L && grid_matches) return(source_path)
  if (volume_count == 1L) {
    ftq_check_overwrite(c(output_path, volume_csv), overwrite)
    values <- as.array(RNifti::readNifti(source_path, internal = FALSE))
    if (!all(is.finite(values))) {
      stop(context, " contains non-finite values.", call. = FALSE)
    }
    # RNifti arrays retain source-header attributes that can override the
    # template affine during writing. Keep only the voxel array and dimensions.
    attributes(values) <- list(dim = dim(values))
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    RNifti::writeNifti(
      values,
      output_path,
      template = write_template,
      datatype = "auto"
    )
    ftq_regional_study_write_image_provenance(
      volume_csv = volume_csv,
      source_path = source_path,
      bvals = bvals,
      volume_index = 1L,
      b_value = NA_real_,
      low_b_max = low_b_max,
      source_grid_matches_region_masks = grid_matches,
      reheadered_to_region_grid = TRUE,
      preparation = "reheader_without_resampling"
    )
    return(output_path)
  }
  if (length(nonspatial_dims) != 1L) {
    stop(
      context, " must be a 3-D image or a 4-D image with one volume dimension.",
      call. = FALSE
    )
  }
  if (is.null(bvals)) {
    stop(
      context, " is four-dimensional; supply matching bvals for this subject.",
      call. = FALSE
    )
  }
  values <- if (is.character(bvals) && length(bvals) == 1L) {
    scan(bvals, quiet = TRUE)
  } else {
    as.numeric(bvals)
  }
  if (length(values) != volume_count || any(!is.finite(values))) {
    stop(
      context, " has ", volume_count,
      " volumes, but its bvals contain ", length(values),
      " finite values.",
      call. = FALSE
    )
  }
  low_b_indices <- which(abs(values) <= low_b_max)
  if (!length(low_b_indices)) {
    stop(
      context, " has no volumes with abs(b-value) <= ", low_b_max, ".",
      call. = FALSE
    )
  }

  ftq_check_overwrite(c(output_path, volume_csv), overwrite)
  low_b <- as.array(RNifti::readNifti(
    source_path,
    volumes = low_b_indices,
    internal = FALSE
  ))
  spatial_dims <- as.integer(source_dims[1:3])
  mean_low_b <- if (length(dim(low_b)) == 3L) {
    low_b
  } else {
    array(
      rowMeans(matrix(low_b, nrow = prod(spatial_dims))),
      dim = spatial_dims
    )
  }
  if (!all(is.finite(mean_low_b))) {
    stop(context, " produced non-finite values in its mean low-b image.",
         call. = FALSE)
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  RNifti::writeNifti(
    mean_low_b,
    output_path,
    template = write_template,
    datatype = "float"
  )
  ftq_regional_study_write_image_provenance(
    volume_csv = volume_csv,
    source_path = source_path,
    bvals = bvals,
    volume_index = low_b_indices,
    b_value = values[low_b_indices],
    low_b_max = low_b_max,
    source_grid_matches_region_masks = grid_matches,
    reheadered_to_region_grid = !grid_matches,
    preparation = "mean_low_b"
  )
  output_path
}

ftq_regional_study_write_image_provenance <- function(
    volume_csv,
    source_path,
    bvals,
    volume_index,
    b_value,
    low_b_max,
    source_grid_matches_region_masks,
    reheadered_to_region_grid,
    preparation) {
  utils::write.csv(
    data.frame(
      source_path = source_path,
      source_md5 = unname(tools::md5sum(source_path)),
      bvals_source = ftq_regional_study_bval_source(bvals),
      bvals_source_md5 = ftq_regional_study_checksum(
        ftq_regional_study_bval_source(bvals)
      ),
      volume_index = volume_index,
      b_value = b_value,
      low_b_max = low_b_max,
      source_grid_matches_region_masks = source_grid_matches_region_masks,
      reheadered_to_region_grid = reheadered_to_region_grid,
      preparation = preparation,
      stringsAsFactors = FALSE
    ),
    volume_csv,
    row.names = FALSE
  )
}

ftq_regional_study_subject_setting <- function(x, subject, subjects) {
  if (is.null(x)) return(NULL)
  nms <- names(x)
  if (!is.null(nms) && all(subjects %in% nms)) return(x[[subject]])
  x
}

ftq_regional_study_scalar <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  as.character(x)[1]
}

ftq_regional_study_bval_source <- function(x) {
  if (is.character(x) && length(x) == 1L) return(x)
  NA_character_
}

ftq_regional_study_checksum <- function(path) {
  if (is.null(path) || !length(path) || is.na(path[1]) ||
      !file.exists(path[1])) {
    return(NA_character_)
  }
  unname(tools::md5sum(path[1]))
}

ftq_regional_study_run_row <- function(subject,
                                       focal_region,
                                       run,
                                       exploratory,
                                       viewer,
                                       started_at,
                                       error,
                                       manifest = NULL,
                                       stage_errors = list(),
                                       requested = NULL) {
  completed <- !is.null(run)
  roi <- if (completed && is.data.frame(run$roi_summary) &&
             nrow(run$roi_summary)) {
    run$roi_summary[1, , drop = FALSE]
  } else {
    NULL
  }
  value <- function(x, default = NA_character_) {
    if (is.null(x) || !length(x) || is.na(x[1])) return(default)
    x[1]
  }
  roi_value <- function(name) {
    if (is.null(roi) || !name %in% names(roi)) return(NA_real_)
    as.numeric(roi[[name]][1])
  }
  provenance_value <- function(name, default) {
    run_manifest <- if (completed) run$manifest else manifest
    manifest_value <- value(
      ftq_value(run_manifest[[name]], run_manifest$analysis[[name]]),
      default
    )
    if (!completed) return(manifest_value)
    if (is.null(roi) || !name %in% names(roi)) return(manifest_value)
    value(roi[[name]], manifest_value)
  }
  output_value <- function(name) {
    if (!completed) return(NA_character_)
    as.character(value(run$outputs[[name]]))
  }
  package_identity <- if (completed && is.list(run$package_identity)) {
    run$package_identity
  } else {
    ftq_package_identity()
  }
  package_value <- function(name) {
    as.character(value(package_identity[[name]]))
  }
  stage_error <- function(stage) {
    as.character(value(stage_errors[[stage]], ""))
  }
  stage_status <- function(stage, object) {
    if (nzchar(stage_error(stage))) return("failed")
    if (!is.null(object)) return("completed")
    if (!is.null(requested) && !isTRUE(requested[[stage]])) return("not_requested")
    if (!completed) return("not_run")
    "not_generated"
  }
  data.frame(
    subject = subject,
    focal_region = focal_region,
    run_status = if (completed && is.null(error)) "completed" else "failed",
    validation_passed = if (completed) isTRUE(run$validation$passed) else FALSE,
    biological_qc_status = if (completed) {
      as.character(value(run$validation$biological_qc_status, "unavailable"))
    } else {
      "unavailable"
    },
    processed_streamlines = roi_value("processed_streamlines"),
    unique_focal_intersecting_streamlines =
      roi_value("unique_focal_intersecting_streamlines"),
    selected_track_space = if (completed) {
      as.character(value(roi$selected_track_space))
    } else {
      NA_character_
    },
    selected_axis_flip = if (completed) {
      as.character(value(roi$selected_axis_flip))
    } else {
      NA_character_
    },
    primary_report = if (completed && !is.null(run$report)) {
      as.character(value(run$report$html))
    } else {
      NA_character_
    },
    exploratory_report = if (!is.null(exploratory)) {
      as.character(value(exploratory$html))
    } else {
      NA_character_
    },
    interactive_3d_qc = if (!is.null(viewer)) {
      as.character(value(viewer$html))
    } else {
      NA_character_
    },
    result_csv = output_value("result_csv"),
    rsc_csv = output_value("rsc_csv"),
    fate_csv = output_value("fate_csv"),
    roi_summary_csv = output_value("roi_summary_csv"),
    endpoint_summary_csv = output_value("endpoint_summary_csv"),
    endpoint_penetration_csv = output_value("endpoint_penetration_csv"),
    endpoint_context_csv = output_value("endpoint_context_csv"),
    local_search_assignments_csv = output_value(
      "local_search_assignments_csv"
    ),
    metadata_json = output_value("metadata_json"),
    validation_json = output_value("validation_json"),
    error = if (is.null(error)) "" else as.character(error),
    started_at = started_at,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    target_intersection_mode = as.character(provenance_value(
      "target_intersection_mode",
      NA_character_
    )),
    local_search_radius_mm = as.numeric(provenance_value(
      "local_search_radius_mm",
      NA_real_
    )),
    package_version = package_value("version"),
    package_remote_type = package_value("remote_type"),
    package_remote_url = package_value("remote_url"),
    package_remote_ref = package_value("remote_ref"),
    package_remote_sha = package_value("remote_sha"),
    package_built = package_value("built"),
    package_loaded_path = package_value("loaded_path"),
    analysis_status = if (completed) "completed" else "failed",
    validation_status = if (!completed) "not_run" else
      if (isTRUE(run$validation$passed)) "passed" else "failed",
    primary_report_status = stage_status("primary_report", run$report),
    exploratory_report_status = stage_status("exploratory_report", exploratory),
    interactive_3d_qc_status = stage_status("interactive_3d_qc", viewer),
    analysis_error = stage_error("analysis"),
    validation_error = stage_error("validation"),
    primary_report_error = stage_error("primary_report"),
    exploratory_report_error = stage_error("exploratory_report"),
    interactive_3d_qc_error = stage_error("interactive_3d_qc"),
    focal_assignment_mode = as.character(provenance_value(
      "focal_assignment_mode", NA_character_
    )),
    stringsAsFactors = FALSE
  )
}

ftq_regional_study_check_labels <- function(labels, context, reserved = character()) {
  safe <- ftq_sanitize_label(labels)
  # A bundle must remain unambiguous on case-insensitive Windows/macOS volumes.
  key <- tolower(sub("[. ]+$", "", safe))
  windows_device <- grepl("^(con|prn|aux|nul|com[1-9]|lpt[1-9])($|[.])", key)
  if (any(!nzchar(key)) || any(safe %in% c(".", "..")) ||
      any(windows_device) || any(grepl("[. ]$", safe)) ||
      any(key %in% tolower(reserved)) || anyDuplicated(key)) {
    stop(context, " labels must produce distinct, non-reserved output folders ",
         "after sanitization (including case-insensitive filesystems and ",
         "Windows device names or trailing periods).", call. = FALSE)
  }
  invisible(safe)
}

ftq_write_regional_study_index <- function(run_index,
                                           path,
                                           output_dir,
                                           description,
                                           example_validation = NULL) {
  esc <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "unavailable"
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    gsub("'", "&#39;", x, fixed = TRUE)
  }
  field <- function(row, name, default = "unavailable") {
    x <- row[[name]]
    if (is.null(x) || !length(x) || is.na(x[1L]) || !nzchar(as.character(x[1L]))) {
      return(default)
    }
    as.character(x[1L])
  }
  link <- function(target, label) {
    if (!length(target) || is.na(target) || !nzchar(target)) return("")
    href <- ftq_report_href(target, dirname(ftq_report_normalize_path(path)))
    paste0("<a href='", esc(href), "'>", esc(label), "</a>")
  }
  mode_text <- function(row) {
    mode <- field(row, "target_intersection_mode", "not recorded")
    labels <- c(full_streamline = "Full streamline", endpoint_voxels = "Endpoint voxels",
                local_search = "Local search")
    text <- if (mode %in% names(labels)) paste0(labels[[mode]], " (", mode, ")") else mode
    if (identical(mode, "local_search")) {
      radius <- suppressWarnings(as.numeric(field(row, "local_search_radius_mm", NA_character_)))
      text <- paste0(text, "; radius ", if (is.finite(radius))
        paste0(format(radius, trim = TRUE), " mm") else "unavailable")
    }
    allocation <- field(row, "focal_assignment_mode", "")
    if (nzchar(allocation)) text <- paste0(text, "; focal allocation: ", allocation)
    text
  }
  modes <- vapply(seq_len(nrow(run_index)), function(i) mode_text(run_index[i, , drop = FALSE]),
                  character(1))
  shared_mode <- length(unique(modes)) == 1L
  download_labels <- c(
    result_csv = "Complete per-voxel data (CSV)",
    voxel_data_csv = "Complete per-voxel data with coordinates (CSV)",
    rsc_csv = "Complete RSC summary (CSV)",
    fate_csv = "Voxel-weighted fate (CSV)",
    roi_summary_csv = "Detailed run summary and unique-streamline fate (CSV)",
    endpoint_summary_csv = "Endpoint assignment audit (CSV)",
    endpoint_penetration_csv = "Target penetration (CSV)",
    endpoint_context_csv = "Endpoint coverage (CSV)",
    local_search_assignments_csv = "Local-search assignments (CSV)",
    metadata_json = "Methods and provenance (JSON)",
    validation_json = "Computational validation (JSON)"
  )
  rows <- vapply(seq_len(nrow(run_index)), function(i) {
    row <- run_index[i, , drop = FALSE]
    legacy_completed <- identical(field(row, "run_status"), "completed")
    analysis <- field(row, "analysis_status", if (legacy_completed) "completed" else "unavailable")
    validation <- field(row, "validation_status", if (isTRUE(row$validation_passed)) "passed" else "not passed")
    primary <- field(row, "primary_report", "")
    action <- if (nzchar(primary)) link(primary, "Results") else paste0(
      "<a href='#run-", i, "'>Results and available files</a>"
    )
    auxiliary <- c(link(field(row, "exploratory_report", ""), "Target report"),
                   link(field(row, "interactive_3d_qc", ""), "3D QC"))
    auxiliary <- auxiliary[nzchar(auxiliary)]
    downloads <- vapply(names(download_labels), function(name) {
      target <- field(row, name, "")
      if (!nzchar(target)) return("")
      # Historical summary CSVs can contain original absolute paths after a
      # move. Do not advertise stale download links; repaired reports remain.
      if (!file.exists(target)) return("")
      paste0("<li>", link(target, download_labels[[name]]), "</li>")
    }, character(1))
    status_fields <- c(primary_report = "Main report", exploratory_report = "Target report",
                       interactive_3d_qc = "3D QC")
    presentation <- vapply(names(status_fields), function(name) {
      status <- field(row, paste0(name, "_status"), "")
      if (!nzchar(status)) return("")
      paste0(esc(status_fields[[name]]), ": ", esc(gsub("_", " ", status, fixed = TRUE)))
    }, character(1))
    presentation <- presentation[nzchar(presentation)]
    error <- field(row, "error", "")
    detail <- paste0(
      "<details id='run-", i, "'", if (nzchar(error) || !nzchar(primary)) " open" else "", ">",
      "<summary>Supporting views, downloads and status</summary>",
      if (length(presentation)) paste0("<p>", paste(presentation, collapse = "; "), ".</p>") else "",
      if (length(auxiliary)) paste0("<p>", paste(auxiliary, collapse = " &middot; "), "</p>") else "",
      if (any(nzchar(downloads))) paste0("<ul>", paste(downloads, collapse = ""), "</ul>") else
        "<p class='muted'>No additional data links are available in this index.</p>",
      "</details>"
    )
    paste0(
      "<tr><th scope='row'>", esc(row$subject), "</th><td>", esc(row$focal_region), "</td>",
      "<td class='", if (analysis == "completed") "pass" else "fail", "'>", esc(analysis), "</td>",
      "<td class='", if (validation == "passed") "pass" else "fail", "'>", esc(validation), "</td>",
      "<td>", esc(gsub("_", " ", field(row, "biological_qc_status"), fixed = TRUE)), "</td>",
      if (!shared_mode) paste0("<td>", esc(modes[[i]]), "</td>") else "",
      "<td>", action, detail, "</td></tr>",
      if (nzchar(error)) paste0("<tr><td colspan='", if (shared_mode) 6L else 7L,
        "' class='fail'>Attention: ", esc(error), "</td></tr>") else ""
    )
  }, character(1))
  title <- if (nrow(run_index) == 1L) {
    paste0("FibrTraQ Results: ", run_index$subject[[1L]], " / ", run_index$focal_region[[1L]])
  } else "FibrTraQ Results"
  guide <- ftq_write_output_guide(output_dir)
  root_files <- c("fibrtraq-run-index.csv" = "Run index (CSV)",
                  "fibrtraq-prepared-inputs.csv" = "Prepared-input provenance (CSV)")
  root_links <- vapply(names(root_files), function(name) {
    target <- file.path(output_dir, name)
    if (file.exists(target)) link(target, root_files[[name]]) else ""
  }, character(1))
  root_links <- c(root_links[nzchar(root_links)],
                  link(guide, "Output definitions and guide"))
  example <- character()
  if (!is.null(example_validation)) {
    counts <- function(x) paste0(sum(x$passed, na.rm = TRUE), "/", nrow(x))
    checked <- isTRUE(example_validation$package_validation_passed) &&
      nrow(example_validation$tables) > 0L && nrow(example_validation$claims) > 0L &&
      isTRUE(all(example_validation$tables$passed)) && isTRUE(all(example_validation$claims$passed))
    example <- c(
      "<section><h2>Bundled example verification</h2>",
      paste0("<p class='", if (checked) "pass" else "fail", "'>Known-answer comparisons: ",
        counts(example_validation$tables), " tables and ", counts(example_validation$claims),
        " claims passed. Computational validation: ",
        if (isTRUE(example_validation$package_validation_passed)) "passed" else "not passed", ".</p>"),
      paste0("<p>Generated HTML files present and nonempty: ", counts(example_validation$reports),
        ". This checks file generation; browser rendering and interactive behavior require separate review.</p>"),
      "<details><summary>Example comparison details</summary>"
    )
    for (family in c("tables", "claims", "reports")) {
      checks <- example_validation[[family]]
      label <- if (family == "claims") checks$claim else checks$artifact
      example <- c(example, paste0("<h3>", esc(family), "</h3><ul>"),
        vapply(seq_len(nrow(checks)), function(i) paste0("<li>", esc(label[[i]]), ": ",
          if (isTRUE(checks$passed[[i]])) "passed" else "FAILED", " \u2014 ",
          esc(checks$details[[i]]), "</li>"), character(1)), "</ul>")
    }
    example <- c(example, "</details></section>")
  }
  html <- c(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    paste0("<title>", esc(title), "</title>"),
    paste0("<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;",
      "max-width:1200px;margin:32px auto;padding:0 20px;color:#1f2933;line-height:1.5;overflow-wrap:anywhere}",
      ".table-scroll{width:100%;overflow-x:auto;-webkit-overflow-scrolling:touch}",
      ".table-scroll:focus{outline:2px solid #2563eb;outline-offset:2px}",
      "table{border-collapse:collapse;width:100%;min-width:700px}",
      "th,td{border:1px solid #d8dee4;padding:10px;text-align:left;vertical-align:top}",
      "thead th{background:#f3f4f6}tbody th{font-weight:600}.pass{color:#176b35}",
      ".fail{color:#9f1239;font-weight:700}.muted{color:#596579}a{color:#175db0}",
      "details{margin-top:10px}summary{cursor:pointer}ul{padding-left:20px}",
      "@media(max-width:700px){body{margin:16px 0;padding:0 12px}h1{font-size:1.6rem}",
      "th,td{padding:6px 8px;font-size:.9rem}}</style></head><body>"),
    paste0("<h1>", esc(title), "</h1>"),
    paste0("<p>", nrow(run_index), if (nrow(run_index) == 1L) " run recorded. " else " runs recorded. ",
      sum(run_index$validation_passed, na.rm = TRUE), " passed computational validation.</p>"),
    if (!nrow(run_index)) "<p>Inputs are being prepared; no analysis run has finished yet.</p>" else "",
    if (shared_mode) paste0("<p><strong>Method:</strong> ", esc(modes[[1L]]), ".</p>") else "",
    "<p>Computational validation checks numerical consistency. Biological QC is a separate review of anatomical context and alignment; unavailable or unreviewed QC is not evidence of anatomical validity.</p>",
    example,
    "<div class='table-scroll' role='region' tabindex='0'>",
    "<table><thead><tr><th>Subject</th><th>Focal region</th><th>Analysis</th><th>Validation</th><th>Biological QC</th>",
    if (!shared_mode) "<th>Method</th>" else "",
    "<th>Results</th></tr></thead><tbody>", rows, "</tbody></table></div>",
    if (length(root_links)) paste0("<h2>Files and definitions</h2><p>",
      paste(root_links, collapse = " &middot; "), "</p>") else "",
    "<p class='muted'>Preparation records retain optional inputs as unavailable when none were supplied. A missing optional reference or b-value file is not an analysis failure. Keep the complete results folder together and extract it before opening the HTML.</p>",
    "</body></html>"
  )
  writeLines(enc2utf8(html), path, useBytes = TRUE)
  invisible(path)
}

ftq_regional_study_relative_path <- function(path, root) {
  ftq_report_relative_path(path, root)
}
