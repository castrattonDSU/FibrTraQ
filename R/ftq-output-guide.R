# Internal output dictionary. Keep definitions tied to the table writers and
# verify coverage against their generated schemas in test_output_guide.R.
ftq_output_field_dictionary <- local({
  # Definitions are fixed for the loaded package. Build once per namespace;
  # ordinary R copy-on-modify keeps callers from changing this private value.
  dictionary <- NULL
  function() {
    if (is.null(dictionary)) dictionary <<- ftq_build_output_field_dictionary()
    dictionary
  }
})

ftq_build_output_field_dictionary <- function() {
  rows <- list()
  add <- function(family, column, definition, units = "text",
                  population = "one run", denominator = "not applicable",
                  missing_zero = "NA means unavailable or not applicable.",
                  method = "selected assignment mode", alias_of = "none") {
    rows[[length(rows) + 1L]] <<- data.frame(
      family = family, column = column, definition = definition, units = units,
      population = population, denominator = denominator,
      missing_zero = missing_zero, method = method, alias_of = alias_of,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  copy <- function(from, to, columns = NULL) {
    source <- do.call(rbind, rows)
    source <- source[source$family == from, , drop = FALSE]
    if (!is.null(columns)) source <- source[source$column %in% columns, , drop = FALSE]
    source$family <- to
    rows[[length(rows) + 1L]] <<- source
  }
  counts <- function(family, definitions, population, method = "selected assignment mode") {
    for (nm in names(definitions)) add(
      family, nm, definitions[[nm]], "count", population,
      missing_zero = "0 means no counted events; NA means unavailable.", method = method
    )
  }
  coords <- function(family, columns = c("x", "y", "z")) {
    for (i in seq_along(columns)) add(
      family, columns[[i]], paste0("One-based focal-image array index along ",
                                  c("x", "y", "z")[[i]],
                                  "; not a scanner-space distance."),
      "voxel index", "one focal voxel", missing_zero = "Positive integer; zero and NA are invalid.",
      method = "focal image grid"
    )
  }
  fraction <- function(family, column, numerator, denominator,
                       population, scale = 100, zero = "NA when denominator is zero.",
                       method = "selected assignment mode") {
    add(family, column, paste0(if (scale == 100) "100 * " else "", numerator,
                              " / ", denominator, "."),
        if (scale == 100) "percent" else "fraction", population, denominator,
        paste("0 means no numerator events when evaluable;", zero), method)
  }

  add("per_voxel", "fNode_Coords", "Coordinate key x_y_z, using one-based focal-image array indices.",
      "voxel index key", "one focal voxel", missing_zero = "Every focal voxel has a key, including zero-assignment voxels.")
  counts("per_voxel", c(
    tracks = "Streamlines allocated to this focal voxel; each streamline contributes at most once to a voxel.",
    local = "Allocated streamlines whose entire selected voxelized path is inside the configured local ROI. Local status can overlap target-positive status.",
    leave_no_target = "Allocated streamlines that are nonlocal and receive no target under the selected assignment rule.",
    `{target}` = "Allocated streamlines assigned to the named target. The target name is supplied by the user; counts across target columns may overlap."
  ), "one focal voxel")
  add("per_voxel", "none", "1 when tracks is zero, otherwise 0. This is zero assignment, not necessarily absence of a traversing streamline under first_intersection.",
      "binary indicator", "one focal voxel", missing_zero = "0 means at least one assigned streamline; 1 means none assigned.")
  add("per_voxel", "cNode_Coords", "Optional legacy coordinate annotation; not a target-count column.",
      "coordinate annotation", "one focal voxel", method = "legacy input only")
  copy("per_voxel", "all_voxels")
  add("all_voxels", "{coordinate}", "Added x/y/z array index. Consult the companion all-voxels-columns.csv: existing target names are retained; conflicting coordinates use voxel_x, voxel_y, voxel_z and a unique suffix.",
      "one-based voxel index", "one focal voxel", missing_zero = "Positive integer; never infer the role from a column name.", method = "focal image grid")
  for (nm in c("column", "role", "source_column", "units", "description")) add(
    "all_voxels_columns", nm, switch(nm,
      column = "Exact column name in the companion all-voxels.csv.",
      role = "Explicit role: coordinate, identifier, count or target. Use this to distinguish coordinates from targets named x, y or z.",
      source_column = "Source field in the original per-voxel table; coordinates are derived from fNode_Coords.",
      units = "Measurement or indexing unit for this column.",
      description = "Human-readable definition of this column."),
    population = "one exported column", method = "additive readable export"
  )
  copy("per_voxel", "local_fraction", "fNode_Coords")
  fraction("local_fraction", "local_fraction", "local", "tracks", "one focal voxel", 1)
  add("local_fraction", "has_assigned_streamlines", "TRUE when tracks > 0; identifies where local_fraction is evaluable.",
      "logical", "one focal voxel", missing_zero = "FALSE means zero assignment and unevaluable fraction.")

  add("rsc", "node", "Target label; the RSC source CSV contains target rows.", population = "one target")
  counts("rsc", c(
    target_associations = "Sum of this target's per-voxel counts; one streamline may contribute at several focal voxels.",
    focal_voxel_streamline_target_associations = "Same count as target_associations.",
    tracks = "Legacy alias of target_associations; not a unique streamline total.",
    voxels_connected = "Focal voxels with a positive count for this target.",
    total_focal_voxels = "All focal voxels, including those with zero assigned streamlines.",
    total_origin_voxels = "Legacy alias of total_focal_voxels."
  ), "one target")
  fraction("rsc", "rsc", "voxels_connected", "total_focal_voxels", "one target", 1,
           "The legacy empty-table case stores 0; it is unevaluable.")
  fraction("rsc", "rsc_percent", "voxels_connected", "total_focal_voxels", "one target", 100,
           "The legacy empty-table case stores 0; it is unevaluable.")
  add("fate", "fate", "Exclusive voxel-weighted class: local, leave_no_target or leave_with_target.", population = "one fate class")
  counts("fate", c(tracks = "Sum of per-voxel counts for this fate; leave_with_target = tracks - local - leave_no_target."),
         "focal-voxel/streamline associations")
  fraction("fate", "percent", "tracks in this row", "sum of tracks across the three fate rows",
           "focal-voxel/streamline associations", zero = "Legacy zero-denominator storage is 0, an unevaluable value.")

  counts("roi_summary", c(
    total_input_streamlines = "Tractogram's declared total, if known; it can differ from the number processed when a limit is used.",
    processed_streamlines = "Streamlines processed during the tractogram pass, before focal-positive filtering.",
    unique_focal_intersecting_streamlines = "Unique streamlines with at least one selected focal voxel.",
    unique_local_streamlines = "Focal-positive streamlines entirely local, irrespective of target assignment.",
    unique_leave_no_target_streamlines = "Focal-positive nonlocal streamlines without a target assignment.",
    unique_leave_with_target_streamlines = "Focal-positive nonlocal streamlines with one or more assigned targets.",
    unique_focal_intersecting_with_target_streamlines = "Focal-positive streamlines with any target assignment, including local ones.",
    unique_local_target_positive_streamlines = "Local and target-positive focal streamlines: one of four exclusive local-by-target groups.",
    unique_local_target_negative_streamlines = "Local and target-negative focal streamlines: one of four exclusive local-by-target groups.",
    unique_nonlocal_target_positive_streamlines = "Nonlocal and target-positive focal streamlines; same as unique_leave_with_target_streamlines.",
    unique_nonlocal_target_negative_streamlines = "Nonlocal and target-negative focal streamlines; same as unique_leave_no_target_streamlines.",
    unique_any_target_positive_streamlines = "Focal-positive streamlines with any assigned target; same as unique_focal_intersecting_with_target_streamlines.",
    unique_any_target_negative_streamlines = "Focal-positive streamlines without any assigned target, combining local and nonlocal groups.",
    unique_exclusive_fate_sum = "Sum of the four exclusive local-by-target groups; should equal unique_focal_intersecting_streamlines.",
    total_focal_voxel_streamline_associations = "Sum of per-voxel tracks, with possible repeat contributions by one streamline across focal voxels.",
    total_focal_voxel_streamline_target_associations = "Sum over all focal voxels and targets; a streamline can contribute across both voxels and targets.",
    `unique_target_streamlines:{target}` = "Unique focal-positive streamlines assigned to this target under the primary method.",
    `unique_full_path_target_streamlines:{target}` = "Unique focal-positive streamlines whose selected full path enters this exact target mask; additive comparator in endpoint/local-search runs."
  ), "one run; unique counts are within that run")
  add("roi_summary", "unique_exclusive_fate_sum_matches_focal", "Arithmetic equality of unique_exclusive_fate_sum and unique_focal_intersecting_streamlines.", "logical")
  settings <- c(
    focal_assignment_mode = "all_focal_voxels counts each occupied focal voxel; first_intersection allocates a streamline only to its first selected focal voxel.",
    streamline_voxelization_mode = "sample_points uses stored tractogram points; segment_interpolation densifies consecutive segments before voxelization.",
    target_intersection_mode = "Primary assignment: full_streamline, endpoint_voxels or local_search. These change target assignment, not the separate full-path endpoint audit.",
    local_roi_defaulted_to_focal = "Whether the focal mask itself defined the local ROI; FALSE means a separately supplied or prepared local region was used.",
    selected_track_space = "Coordinate mapping actually used for tractogram points.",
    track_voxel_size_mm = "Comma-separated TrackVis voxel sizes used for coordinate conversion, in millimetres.",
    selected_axis_flip = "Comma-separated logical axis flips in x,y,z order.",
    selected_discretization = "Ordered voxel discretization candidates. For any, a focal-positive path uses the first candidate that qualifies it.",
    reverse_streamlines = "Whether stored point order was reversed before counting.",
    backend = "Counting implementation used: R or compiled backend; not a biological method.",
    mapping_selection_rule = "Rule selecting automatic coordinate mapping; target/local hits do not select or break ties.",
    local_search_spatial_units = "Spatial-unit declaration in the reference image header.",
    local_search_spatial_units_interpretation = "Recorded interpretation of those units, including an explicit millimetre assumption for unknown units.",
    summary_available = "FALSE in a fallback record when a reloaded per-voxel CSV lacks unique-streamline provenance.",
    note = "Explanation of unavailable summary provenance."
  )
  for (nm in names(settings)) add("roi_summary", nm, settings[[nm]])
  add("roi_summary", "local_search_radius_mm", "Inclusive maximum endpoint-to-labelled-voxel-centre distance.", "mm", method = "local_search only")
  add("roi_summary", "local_search_spatial_units_scale_to_mm", "Multiplicative factor converting declared image-space distance to millimetres.", "conversion factor", method = "local_search only")

  add("connection_summary", "cNode", "Target label, or local for the configured local region. Local can overlap target-positive rows.", population = "one target/local row")
  counts("connection_summary", c(
    unique_streamline_connections = "ROI-level unique streamline count for this target or local class; NA if unavailable from saved provenance.",
    connected_fNode_voxels = "Focal voxels with a positive target/local count.",
    total_fNode_voxels = "All focal voxels, including zero-assignment voxels."
  ), "one target/local row")
  fraction("connection_summary", "rsc_percent", "connected_fNode_voxels", "total_fNode_voxels", "one target/local row")
  add("report_fate", "fate", "Exclusive class: local, leave_no_target or leave_with_target. Read summary_basis before interpreting counts.", population = "one fate class")
  counts("report_fate", c(count = "Count for this class in the stated summary_basis.",
                          denominator_count = "All focal-positive unique streamlines, or all focal-voxel/streamline associations in the explicit fallback."), "summary_basis")
  add("report_fate", "summary_basis", "unique_focal_intersecting_streamlines when available; focal_voxel_streamline_associations is a clearly declared fallback.")
  fraction("report_fate", "percent_of_denominator", "count", "denominator_count", "summary_basis",
           zero = "Legacy zero-denominator storage is 0; the fraction is unevaluable.")
  copy("per_voxel", "top_voxels", c("tracks", "local", "leave_no_target", "{target}"))
  coords("top_voxels")
  counts("top_voxels", c(leave_with_target = "tracks - local - leave_no_target; nonlocal allocated streamlines with any assigned target."), "one focal voxel")
  counts("top_voxels", c(`{derived_leave_with_target}` = "Same derived nonlocal target-positive count as tracks - local - leave_no_target. If a genuine target is named leave_with_target, the derived field uses derived_leave_with_target, repeatedly adding derived_ until unique; preserve the target's original label."), "one focal voxel")
  add("top_voxels", "leave_no_target.1", "Duplicate leave_no_target column in older top-50 previews; retained in historical files, omitted from newly generated ranked exports.", "count", "one focal voxel", alias_of = "leave_no_target", method = "historical preview only")

  endpoint_method <- "additive full-path audit; exact or Chebyshev one-voxel dilation, independent of local_search radius"
  add("endpoint_summary", "target", "Target whose exact full-path intersection population is audited.", population = "one target/tolerance", method = endpoint_method)
  add("endpoint_summary", "tolerance", "exact or dilated1; the latter includes the 26 neighbouring voxels and is not a uniform millimetre radius.", method = endpoint_method)
  add("endpoint_summary", "tolerance_definition", "Human-readable definition of the endpoint-membership tolerance.", method = endpoint_method)
  counts("endpoint_summary", c(
    unique_focal_intersecting_streamlines = "All unique focal-positive streamlines in the run.",
    target_intersection_streamlines = "Unique focal-positive streamlines whose selected full path occupies at least one exact target voxel.",
    focal_endpoint_streamlines = "Focal-positive streamlines with at least one focal endpoint at this tolerance, irrespective of target intersection.",
    paired_endpoint_streamlines = "Focal-positive streamlines with opposite focal and target endpoints, including ambiguous shared-mask membership and pairs without an exact target intersection.",
    target_intersection_with_paired_endpoints = "Exact target-intersection streamlines with opposite focal and target endpoints; this broad count can include shared-mask ambiguity.",
    target_intersection_with_any_target_endpoint = "Exact target-intersection streamlines with at least one target endpoint, including ambiguous cases.",
    target_intersection_with_any_focal_endpoint = "Exact target-intersection streamlines with any focal endpoint: shared-mask-ambiguous, paired and focal-only classes combined.",
    focal_endpoint_with_target_intersection_streamlines = "Focal-positive streamlines with a focal endpoint and an exact target intersection, including ambiguous membership.",
    shared_mask_ambiguous = "Target-intersection streamlines with an endpoint in both focal and target masks at this tolerance; precedence class in the five-class partition.",
    paired_focal_target_endpoints = "Unambiguous opposite focal/target endpoints among exact target-intersection streamlines; excludes shared-mask ambiguity.",
    target_endpoint_only = "Unambiguous target endpoint present, no focal endpoint, among target-intersection streamlines.",
    focal_endpoint_only = "Unambiguous focal endpoint present, no target endpoint, among target-intersection streamlines.",
    neither_relevant_endpoint = "Neither relevant endpoint in the five-class target-intersection partition.",
    focal_mask_voxels = "Focal-mask voxel count at the displayed exact/dilated tolerance.",
    target_mask_voxels = "Target-mask voxel count at the displayed exact/dilated tolerance.",
    focal_target_overlap_voxels = "Voxel overlap of focal and target masks at this tolerance."
  ), "one target/tolerance", endpoint_method)
  endpoint_ratios <- list(
    target_intersection_percent_of_unique_focal_intersecting_streamlines = c("target_intersection_streamlines", "unique_focal_intersecting_streamlines"),
    paired_percent_of_target_intersections = c("target_intersection_with_paired_endpoints", "target_intersection_streamlines"),
    paired_percent_of_focal_endpoint_streamlines = c("paired_endpoint_streamlines", "focal_endpoint_streamlines"),
    target_endpoint_percent_of_target_intersections = c("target_intersection_with_any_target_endpoint", "target_intersection_streamlines"),
    focal_endpoint_percent_of_target_intersections = c("target_intersection_with_any_focal_endpoint", "target_intersection_streamlines"),
    target_intersection_percent_of_focal_endpoint_streamlines = c("focal_endpoint_with_target_intersection_streamlines", "focal_endpoint_streamlines"),
    focal_target_overlap_percent_of_smaller_mask = c("focal_target_overlap_voxels", "min(focal_mask_voxels, target_mask_voxels)")
  )
  for (cl in c("shared_mask_ambiguous", "target_endpoint_only", "focal_endpoint_only", "neither_relevant_endpoint")) {
    endpoint_ratios[[paste0(cl, "_percent_of_target_intersections")]] <- c(cl, "target_intersection_streamlines")
  }
  for (nm in names(endpoint_ratios)) fraction("endpoint_summary", nm, endpoint_ratios[[nm]][1],
                                             endpoint_ratios[[nm]][2], "one target/tolerance", method = endpoint_method)
  for (axis in c("x", "y", "z")) add("endpoint_summary", paste0("voxel_size_", axis, "_mm"),
                                        paste("Reference voxel spacing along", axis, "after conversion of declared header units to millimetres, or from an explicit stored millimetre spacing attribute."), "mm", method = endpoint_method)
  copy("roi_summary", "endpoint_summary", c("streamline_voxelization_mode", "target_intersection_mode"))

  copy("endpoint_summary", "endpoint_penetration", c("target", "target_intersection_streamlines", "streamline_voxelization_mode"))
  counts("endpoint_penetration", c(
    unique_target_voxels = "Exact distinct target voxels visited by a selected full path; histogram bin value, not depth in mm.",
    streamline_count = "Unique target-intersection streamlines in this penetration bin."
  ), "one target/penetration bin", endpoint_method)
  fraction("endpoint_penetration", "percent_of_target_intersections", "streamline_count", "target_intersection_streamlines", "one target/penetration bin", method = endpoint_method)
  add("endpoint_penetration", "penetration_interpretation", "Stored-point or segment-interpolated penetration definition; sparse stored points may miss intervening target voxels.", method = endpoint_method)
  copy("endpoint_summary", "penetration_summary", c("target", "target_intersection_streamlines"))
  bins <- c(one_voxel_touch = "exactly 1", two_voxel_only = "exactly 2", three_to_four_voxels = "3 or 4",
            five_to_nine_voxels = "5 through 9", ten_or_more_voxels = "at least 10", at_most_two_voxels = "at most 2",
            at_least_three_voxels = "at least 3", at_least_five_voxels = "at least 5", at_least_ten_voxels = "at least 10")
  for (nm in names(bins)) {
    field <- paste0(nm, "_streamlines")
    counts("penetration_summary", stats::setNames(paste("Target-intersection streamlines visiting", bins[[nm]], "distinct exact target voxels."), field), "one target", endpoint_method)
    fraction("penetration_summary", paste0(nm, "_percent_of_target_intersections"), field, "target_intersection_streamlines", "one target", method = endpoint_method)
  }
  for (q in c("p25", "median", "p75", "p95")) add("penetration_summary", paste0("penetration_", q, "_unique_target_voxels"),
      paste(q, "of unique target voxels, using the first histogram bin whose cumulative streamline count reaches the specified percentile."),
      "voxels", "target-intersection streamlines", missing_zero = "NA when no target-intersection streamlines.", method = endpoint_method)

  context_method <- "all valid mapped streamlines before focal filtering; first discretization candidate"
  add("endpoint_context", "context", "reference_grid or a user-named trusted anatomical context mask.", method = context_method)
  add("endpoint_context", "context_type", "reference_grid checks bounds; binary_mask checks membership in supplied anatomy. These do not establish biological validity.", method = context_method)
  counts("endpoint_context", c(
    valid_streamlines = "Streamlines with valid mapped points, before focal-intersection filtering.",
    endpoint_count = "Available endpoint observations: two per mapped streamline with at least two stored points, one for a single-point record.",
    endpoints_inside = "Endpoints inside this grid/mask.",
    endpoints_outside = "Endpoints not inside this grid/mask; endpoint_count - endpoints_inside.",
    endpoints_nonfinite = "Endpoints whose mapped coordinates are nonfinite (grid row only).",
    endpoints_one_voxel_outside = "Finite endpoints outside bounds with maximum axis excess at most one voxel (grid row only).",
    endpoints_more_than_one_voxel_outside = "Finite endpoints outside bounds with maximum axis excess greater than one voxel (grid row only).",
    streamlines_with_any_endpoint_inside = "Valid mapped streamlines with at least one endpoint inside this grid/mask.",
    streamlines_with_both_endpoints_inside = "Mapped streamlines with all available endpoints inside this grid/mask; a single-point record has one available endpoint.",
    mask_voxels = "Grid volume in voxels or positive anatomical-mask voxels."
  ), "all valid mapped streamlines/endpoints, not just focal-positive ones", context_method)
  fraction("endpoint_context", "endpoints_inside_percent", "endpoints_inside", "endpoint_count", "all valid mapped endpoints", method = context_method)
  fraction("endpoint_context", "streamlines_with_both_endpoints_inside_percent", "streamlines_with_both_endpoints_inside", "valid_streamlines", "all valid mapped streamlines", method = context_method)
  for (axis in c("x", "y", "z")) for (side in c("below", "above")) {
    for (distance in c("grid", "one_voxel", "more_than_one_voxel")) {
      nm <- if (distance == "grid") paste(axis, side, "grid", sep = "_") else paste(axis, distance, side, sep = "_")
      add("endpoint_context", nm, paste("Endpoints", side, "the", axis, "grid boundary; distance category:", gsub("_", " ", distance), ". Axis counts overlap when more than one boundary is crossed."),
          "endpoint count", "all valid mapped endpoints", missing_zero = "0 means none; NA for anatomical-mask rows, where boundary distances are not applicable.", method = context_method)
    }
  }
  for (nm in c("source_path", "checksum_md5", "dimensions", "affine_match_reference")) add("endpoint_context", nm,
      switch(nm, source_path = "Source context mask path; NA for the reference-grid row.", checksum_md5 = "MD5 checksum of the source context mask; NA for the reference-grid row.",
             dimensions = "Grid dimensions in x by y by z array order.", affine_match_reference = "Whether the context grid matches the focal/reference affine; the reference-grid row is TRUE."), method = context_method)

  add("local_search_assignments", "streamline_index", "One-based record position in the input tractogram, not a rank or biological subject identifier.", "streamline index", "one focal-positive streamline", method = "local_search")
  for (end in c("start", "end")) {
    add("local_search_assignments", paste0(end, "_status"), "focal, a target label, no_node_within_radius, ambiguous_nearest_node or endpoint_unavailable. Stored point order defines start/end.", population = "one focal-positive streamline endpoint", method = "local_search")
    add("local_search_assignments", paste0(end, "_distance_mm"), "Nearest-node scanner-space Euclidean distance from the precise endpoint to a labelled voxel centre; ambiguity can retain a finite distance.", "mm", "one focal-positive streamline endpoint", missing_zero = "0 means endpoint coincides with a labelled voxel centre; NA when no eligible node/distance is available.", method = "local_search")
    add("local_search_assignments", paste0(end, "_ambiguous"), "Whether competing distinct labels tied or overlapped at the nearest candidate.", "logical", "one focal-positive streamline endpoint", method = "local_search")
  }
  add("local_search_assignments", "paired_target", "Assigned target label only when the two endpoints form an unambiguous focal/target pair.", population = "one focal-positive streamline", missing_zero = "NA means no assigned pair; consult endpoint statuses.", method = "local_search")
  for (part in c("focal", "target", "max")) add("local_search_assignments", paste0("paired_", part, "_distance_mm"),
      switch(part, focal = "Distance for the focal-assigned end of a successful pair.", target = "Distance for the target-assigned end of a successful pair.", max = "Maximum of the two successful-pair endpoint distances."),
      "mm", "one focal-positive streamline", missing_zero = "NA if unpaired; 0 is an exact centre match.", method = "local_search")
  add("local_search_assignments", "paired_target_intersects_full_path", "Whether the selected voxelized path also occupies an exact voxel of the assigned target. FALSE for an assigned pair means proximity only.",
      "logical", "one focal-positive streamline", missing_zero = "NA if no pair is assigned.", method = "local_search")
  add("local_search_assignments", "search_radius_mm", "Inclusive configured endpoint search radius, repeated for self-contained rows.", "mm", method = "local_search")
  add("local_search_summary", "target", "Target receiving at least one successful local-search pair; targets with no pairs have no row.", population = "one assigned target", method = "local_search")
  counts("local_search_summary", c(unique_connections = "Successful focal/target pairs; one row per focal-positive streamline in the source assignment table.",
      connections_with_full_path_intersection = "Assigned pairs whose selected full path enters the exact assigned target mask.",
      connections_without_full_path_intersection = "Assigned pairs whose selected full path does not enter the assigned target mask (proximity only)."), "one assigned target", "local_search")
  fraction("local_search_summary", "full_path_intersection_percent_of_connections", "connections_with_full_path_intersection", "unique_connections", "one assigned target", method = "local_search")
  distance_fields <- c(focal_distance_median_mm = "median paired_focal_distance_mm", target_distance_median_mm = "median paired_target_distance_mm",
      paired_max_distance_q25_mm = "25th percentile paired_max_distance_mm", paired_max_distance_median_mm = "median paired_max_distance_mm",
      paired_max_distance_q75_mm = "75th percentile paired_max_distance_mm", paired_max_distance_q95_mm = "95th percentile paired_max_distance_mm",
      paired_max_distance_max_mm = "maximum paired_max_distance_mm",
      full_path_target_distance_median_mm = "median paired_target_distance_mm among pairs with full-path intersection",
      full_path_max_distance_median_mm = "median paired_max_distance_mm among pairs with full-path intersection",
      full_path_max_distance_q95_mm = "95th percentile paired_max_distance_mm among pairs with full-path intersection",
      proximity_only_target_distance_median_mm = "median paired_target_distance_mm among proximity-only pairs",
      proximity_only_max_distance_median_mm = "median paired_max_distance_mm among proximity-only pairs",
      proximity_only_max_distance_q95_mm = "95th percentile paired_max_distance_mm among proximity-only pairs")
  for (nm in names(distance_fields)) add("local_search_summary", nm, paste0(distance_fields[[nm]], "; finite distances only, R quantile type 7."), "mm", "successful assigned pairs in the stated target/subgroup", missing_zero = "NA for an empty subgroup; 0 means exact centre match.", method = "local_search")
  copy("local_search_assignments", "local_search_summary", "search_radius_mm")
  add("local_search_endpoint_status", "endpoint", "start or end in stored streamline point order.", method = "local_search")
  add("local_search_endpoint_status", "status", "Endpoint status label as defined in local_search_assignments.", method = "local_search")
  counts("local_search_endpoint_status", c(streamline_count = "Focal-positive streamlines with this endpoint status.", denominator_focal_intersecting_streamlines = "All rows in the focal-positive local-search assignment table."), "one endpoint/status category", "local_search")
  fraction("local_search_endpoint_status", "percent_of_focal_intersecting_streamlines", "streamline_count", "denominator_focal_intersecting_streamlines", "one endpoint/status category", method = "local_search")

  copy("per_voxel", "exploratory_per_voxel", c("fNode_Coords", "tracks", "local", "leave_no_target"))
  coords("exploratory_per_voxel")
  counts("exploratory_per_voxel", c(leave_with_target = "tracks - local - leave_no_target.", total_target_hits = "Sum of target associations across this voxel's targets; can exceed tracks.",
      target_diversity = "Number of targets with a positive count at this focal voxel.", primary_target_hits = "Largest target association count at this focal voxel."), "one focal voxel")
  add("exploratory_per_voxel", "primary_target", "Label with the largest count; ties choose the first target column. This is an association-based label, not proof of anatomical destination.", population = "one focal voxel", missing_zero = "none when no target associations.")
  add("exploratory_per_voxel", "primary_target_index", "One-based category code into the primary-target LUT; 0 for no target associations. Do not interpret as a continuous measure.", "category index", "one focal voxel", missing_zero = "0 means no target; positive integers refer to LUT index.")
  fraction("exploratory_per_voxel", "primary_target_fraction", "primary_target_hits", "total_target_hits", "one focal voxel", 1, "Stored 0 for no target associations, an unevaluable fraction.")
  fraction("exploratory_per_voxel", "target_hits_per_assigned_track", "total_target_hits", "tracks", "one focal voxel", 1, "Stored 0 for zero tracks; may exceed 1 with multiple targets.")
  add("exploratory_per_voxel", "target_entropy_bits", "-sum(p * log2(p)) over positive target proportions p = target_count / total_target_hits.", "bits", "one focal voxel", "total_target_hits", "0 for a single positive target and also stored 0 for no target associations; use coverage to distinguish them.")
  add("exploratory_per_voxel", "target_entropy_fraction", "target_entropy_bits / log2(target_diversity) when diversity > 1.", "fraction", "one focal voxel", "log2(number of positive targets)", "Stored 0 for diversity 0 or 1; absence and concentration are distinct.")
  for (nm in c("leave_no_target", "leave_with_target")) fraction("exploratory_per_voxel", paste0(nm, "_fraction"), nm, "tracks", "one focal voxel", 1, "Stored 0 when tracks is zero; unevaluable, not evidence of absence.")
  add("exploratory_target_totals", "target", "Target label; sorted by decreasing association total, then target label.", population = "one target")
  counts("exploratory_target_totals", c(focal_voxel_streamline_target_associations = "Sum of this target's per-voxel associations.", focal_voxels_with_association = "Focal voxels with a positive count for this target."), "one target")
  fraction("exploratory_target_totals", "percent_of_focal_voxel_streamline_associations", "focal_voxel_streamline_target_associations", "sum of per-voxel tracks", "one target", zero = "Stored 0 for zero total tracks; percentages across targets need not sum to 100.")
  add("exploratory_lut", "index", "One-based primary_target_index category code; zero is reserved for no target and has no LUT row.", "category index", "one target")
  add("exploratory_lut", "target", "Original target label associated with the category code.", population = "one target")
  add("exploratory_lut", "color", "Display colour for the target category; a visual key, not a measurement.", "colour code", "one target")
  copy("roi_summary", "exploratory_provenance", c("target_intersection_mode", "local_search_radius_mm"))

  # Supporting file-level provenance, preparation and diagnostics.
  text_fields <- c(subject = "Explicit user-supplied subject identifier; retain leading zeros and suffixes.",
      focal_region = "Focal-region label for this run.", session = "Optional session identifier.", description = "User-supplied run description.",
      prefix = "Filename prefix preserving run identity.", run_status = "Analysis completion state; inspect error and separate report status when present.",
      biological_qc_status = "Separate anatomical/endpoint QC interpretation; computational success is not biological validation.",
      error = "Recorded analysis or presentation error, or empty text if none; inspect individual stage errors when available.", started_at = "Recorded start timestamp with timezone offset.", completed_at = "Recorded completion timestamp with timezone offset.")
  for (nm in names(text_fields)) add("run_index", nm, text_fields[[nm]])
  add("run_index", "validation_passed", "Whether executable computational validation passed for the completed run.", "logical")
  copy("roi_summary", "run_index", c("processed_streamlines", "unique_focal_intersecting_streamlines", "selected_track_space", "selected_axis_flip", "target_intersection_mode", "local_search_radius_mm"))
  output_paths <- c(primary_report = "main Results HTML", exploratory_report = "supporting Target report HTML", interactive_3d_qc = "interactive sampled 3-D QC HTML",
      result_csv = "original complete per-voxel table", rsc_csv = "legacy full RSC table", fate_csv = "voxel-weighted fate table", roi_summary_csv = "wide unique-streamline run record",
      endpoint_summary_csv = "additive full-path endpoint audit", endpoint_penetration_csv = "target-penetration histogram", endpoint_context_csv = "endpoint coverage table",
      local_search_assignments_csv = "one-row-per-focal-streamline local-search audit", metadata_json = "machine-readable run provenance", validation_json = "machine-readable validation record")
  for (nm in names(output_paths)) add("run_index", nm, paste("Saved path to", output_paths[[nm]], "; paths document the generating workstation; use portable HTML navigation after moving folders."), "file path", missing_zero = "NA means not generated or unavailable; consult analysis/report status.")
  for (stage in c("analysis", "validation", "primary_report", "exploratory_report", "interactive_3d_qc")) {
    add("run_index", paste0(stage, "_status"), if (stage == "validation") {
      "Computational validation stage: passed, failed or not_run. Independent of report-generation status."
    } else {
      paste0(stage, " stage: completed, failed, not_requested, not_run or not_generated, as applicable. Completed analysis outputs remain available after a presentation failure.")
    })
    add("run_index", paste0(stage, "_error"), paste("Error recorded for the", stage, "stage; empty text if absent."))
  }
  copy("roi_summary", "run_index", "focal_assignment_mode")
  for (nm in c("version", "remote_type", "remote_url", "remote_ref", "remote_sha", "built", "loaded_path")) add("run_index", paste0("package_", nm),
      paste("Loaded package identity:", gsub("_", " ", nm), ". remote_sha identifies source when present; a development version alone does not."))
  copy("run_index", "prepared_inputs", "subject")
  prepared_paths <- c(tractogram = "original tractogram", reference_source = "source reference image", bvals_source = "source b-value file (NA when numeric b-values were supplied in memory)", reference_image = "used/prepared three-dimensional reference", background_source = "source background image", background_image = "used/prepared background", local_mask = "used local-region mask")
  for (nm in names(prepared_paths)) add("prepared_inputs", nm, paste("Path to", prepared_paths[[nm]], "; unavailable optional inputs do not imply failure."), "file path")
  for (nm in c("reference_source_md5", "bvals_source_md5", "background_source_md5")) add("prepared_inputs", nm, "MD5 checksum of the corresponding source file, when supplied.", "checksum")
  add("prepared_inputs", "reheadered_to_region_grid", "Whether explicit, user-authorized grid-header replacement was requested for prepared images; no interpolation is implied.", "logical")
  for (nm in c("source_path", "source_md5", "bvals_source", "bvals_source_md5")) add("prepared_volumes", nm,
      if (grepl("md5", nm)) "MD5 checksum of the corresponding source file; NA when no file was supplied." else "Original image/b-value source path; bvals_source is NA for numeric b-values supplied in memory.", if (grepl("md5", nm)) "checksum" else "file path", "one selected source volume")
  add("prepared_volumes", "volume_index", "One-based source fourth-dimension volume index included in the prepared image.", "volume index", "one selected source volume")
  add("prepared_volumes", "b_value", "Diffusion b-value of the selected source volume.", "s/mm^2", "one selected source volume")
  add("prepared_volumes", "low_b_max", "Inclusive maximum b-value for low-b averaging.", "s/mm^2")
  for (nm in c("source_grid_matches_region_masks", "reheadered_to_region_grid")) add("prepared_volumes", nm,
      if (nm == "source_grid_matches_region_masks") "Whether the source grid already matched regional masks." else "Whether a prepared derivative received the trusted regional-grid header.", "logical")
  add("prepared_volumes", "preparation", "Preparation operation: mean_low_b or the recorded three-dimensional image preparation route.")
  for (nm in c("label", "path", "md5")) add("union_inputs", nm, switch(nm, label = "Regional mask label.", path = "Original regional mask path.", md5 = "MD5 checksum of that source mask."), population = "one union input mask")
  counts("union_inputs", c(nonzero_voxels = "Foreground voxels in this source mask."), "one union input mask")
  for (nm in c("mask_a", "mask_b")) add("union_overlaps", nm, "Label of one member of the mask pair.", population = "one mask pair")
  counts("union_overlaps", c(overlap_voxels = "Foreground voxels shared by both aligned masks; pairs are not independent samples."), "one mask pair")

  counts("tractogram_context", c(processed_streamlines = "All streamlines visited in the context pass, before focal filtering.", declared_streamlines = "Tractogram header total; NA if unknown.",
      nonzero_tractogram_occupied_voxels = "Reference-grid voxels with at least one sampled or interpolated streamline occupancy.", total_streamline_voxel_occupancy_events = "Sum of density; one streamline contributes at most once per voxel.",
      max_density = "Largest streamline occupancy count at any voxel.", total_sample_points = "Number of stored points successfully mapped.", total_voxelization_points = "Number of points after the selected voxelization, including interpolated points when enabled.",
      out_of_bounds_sample_points = "Out-of-bounds voxelization points; despite the legacy name this includes interpolated points in segment_interpolation mode.", mapping_error_streamlines = "Streamlines whose coordinates failed conversion in the context pass."), "whole tractogram context pass", "context density discretization; independent of focal/target assignment")
  add("tractogram_context", "tractogram_occupancy_volume_mm3", "Occupied voxel count multiplied by physical voxel volume. This is tractogram occupancy, not anatomical tissue volume.", "mm^3")
  add("tractogram_context", "note", "Interpretation reminder that density context is not anatomical brain volume.")
  for (space in c("voxel", "world")) for (axis in c("x", "y", "z")) for (edge in c("min", "max")) add("tractogram_context", paste(space, "bbox", axis, edge, sep = "_"),
      paste(edge, axis, "bound of occupied voxels; world bounds use transformed corners of the voxel bounding box."), if (space == "voxel") "one-based voxel index" else "mm", missing_zero = "NA when no in-grid occupied voxel.")

  for (nm in c("mask", "template_image", "reference_mask")) add("mask_grid", nm, "Source mask/template/reference image path for this grid audit.", "file path")
  for (nm in c("mask_dims", "template_dims")) add("mask_grid", nm, "Image array dimensions in x by y by z order.", "voxel dimensions")
  for (nm in c("mask_voxel_size", "template_voxel_size")) add("mask_grid", nm, "Raw NIfTI pixdim values in x,y,z order; interpret using header units.", "header spatial units")
  for (nm in c("dims_match", "voxel_size_match", "affine_match", "grid_match", "reference_grid_match_template")) add("mask_grid", nm,
      switch(nm, dims_match = "Equality of mask and template dimensions.", voxel_size_match = "Equality of raw pixdim values to tolerance 1e-6.", affine_match = "Equality of unit-normalized declared affines.", grid_match = "Both dimensions and affine match.", reference_grid_match_template = "Reference mask shares the template grid."), "logical")
  counts("mask_grid", c(nonzero_voxels = "Finite nonzero voxels in the audited mask.", overlap_voxels = "Mask/reference foreground intersection, computed only on matching grids.", union_voxels = "Mask/reference foreground union, computed only on matching grids."), "one mask/reference comparison")
  add("mask_grid", "reference_label", "Supplied reference-mask label; NA if no reference was supplied.")
  fraction("mask_grid", "jaccard", "overlap_voxels", "union_voxels", "one mask/reference comparison", 1, "NA if absent references, grids differ, or union is empty.")
  for (space in c("voxel", "world")) for (axis in c("x", "y", "z")) for (edge in c("min", "max")) add("mask_grid", paste(space, "bbox", axis, edge, sep = "_"),
      paste(edge, axis, "bound of finite nonzero mask voxels, transformed to world coordinates where indicated."), if (space == "voxel") "one-based voxel index" else "NIfTI world units", missing_zero = "NA for an empty mask.")

  copy("run_index", "validation_suite", c("subject", "session", "description", "prefix", "run_status", "validation_passed", "biological_qc_status", "error", "started_at", "completed_at", "target_intersection_mode", "local_search_radius_mm", names(output_paths)[!names(output_paths) %in% c("primary_report", "exploratory_report", "interactive_3d_qc", "local_search_assignments_csv")]))
  for (nm in c("manifest_path", "dataset_root", "output_dir", "report_html")) add("validation_suite", nm, "Manifest, dataset, output or report location recorded for the validation run.", "file path")
  add("validation_suite", "preflight_passed", "Whether input preflight passed; unavailable if not reached.", "logical")
  counts("validation_suite", c(preflight_errors = "Failed preflight checks of severity error.", preflight_warnings = "Failed preflight checks of severity warning.", preflight_recommendations = "Failed preflight checks of severity recommendation.", focal_voxels = "Complete focal-voxel count.",
      focal_voxel_streamline_associations = "Sum of per-voxel tracks.", local_associations = "Sum of per-voxel local.", leave_no_target_associations = "Sum of per-voxel leave_no_target.", leave_with_target_associations = "Sum of tracks - local - leave_no_target."), "one validation manifest")
  add("validation_suite", "validation_failed_checks", "Semicolon-separated failed computational check names; empty when none.")
  add("validation_suite", "endpoint_analysis", "Whether exact_and_dilated1 or none was selected for the additive audit.")
  add("validation_suite", "focal_voxel_streamline_target_associations", "Semicolon-separated target=count association totals; this field is text rather than one numeric total.")

  # Statistical helper tables are also saved by the installed dataset recipe.
  for (family in c("stat_roi", "stat_targets", "stat_voxels")) {
    for (nm in c("subject", "group", "roi")) add(family, nm,
      switch(nm, subject = "User-supplied subject identifier; preserve leading zeros.", group = "Optional user-supplied group/condition label.", roi = "Optional focal-region label."),
      population = "one subject/focal region")
  }
  counts("stat_roi", c(focal_voxels = "All focal voxels.", connected_focal_voxels = "Focal voxels with tracks > 0, not just target-positive voxels.",
      focal_voxel_streamline_associations = "Sum of per-voxel tracks.", local_associations = "Sum of per-voxel local.", leave_no_target_associations = "Sum of per-voxel leave_no_target.",
      leave_with_target_associations = "Sum of tracks - local - leave_no_target.", focal_voxel_streamline_target_associations = "Sum of all target columns over all focal voxels."), "one subject/focal region; voxel-weighted associations")
  fraction("stat_roi", "leave_no_target_fraction", "leave_no_target_associations", "focal_voxel_streamline_associations", "one subject/focal region", 1)
  fraction("stat_roi", "target_hits_per_focal_voxel_streamline_association", "focal_voxel_streamline_target_associations", "focal_voxel_streamline_associations", "one subject/focal region", 1)
  add("stat_roi", "median_target_diversity", "Median number of positive targets over focal voxels with tracks > 0 (includes target-negative assigned voxels).", "targets", "positive-assignment focal voxels", missing_zero = "NA if no assigned voxel; zero means assigned voxels lack targets.")
  add("stat_roi", "median_entropy_fraction", "Median normalized target entropy over focal voxels with tracks > 0; stored zero for no positive target or only one target is included.", "fraction", "positive-assignment focal voxels", missing_zero = "NA if no assigned voxel; zero can mean no targets or concentrated association.")
  add("stat_targets", "target", "Target label.", population = "one subject/focal/target row")
  counts("stat_targets", c(focal_voxel_streamline_target_associations = "Sum of this target's per-voxel counts.", connected_focal_voxels = "Focal voxels with a positive count for this target."), "one subject/focal/target row")
  fraction("stat_targets", "target_fraction_of_focal_voxel_streamline_associations", "focal_voxel_streamline_target_associations", "sum of per-voxel tracks", "one subject/focal/target row", 1)
  fraction("stat_targets", "target_fraction_of_target_hits", "focal_voxel_streamline_target_associations", "sum of all target associations", "one subject/focal/target row", 1)
  copy("per_voxel", "stat_voxels", c("fNode_Coords", "tracks"))
  coords("stat_voxels")
  copy("exploratory_per_voxel", "stat_voxels", "target_diversity")
  fraction("stat_voxels", "leave_no_target_fraction", "leave_no_target", "tracks", "one focal voxel", 1)
  counts("stat_voxels", c(target_hit_events = "Sum of all target associations at this voxel; an alias of exploratory total_target_hits."), "one focal voxel")
  add("stat_voxels", "entropy_fraction", "Normalized target entropy as in exploratory target_entropy_fraction; stored zero for no target or one target.", "fraction", "one focal voxel", "log2(number of positive targets)")

  for (family in c("target_unique", "target_associations")) {
    add(family, "target", "Target label.", population = "one target")
    add(family, "target_family", "Optional user-supplied target family; defaults to unspecified.", population = "one target")
    counts(family, c(target_mask_voxels = "Optional supplied target-mask foreground size; NA when not provided."), "one target")
  }
  counts("target_unique", c(unique_target_positive_focal_intersecting_streamlines = "ROI-level unique focal-positive streamlines assigned to this target; read from the source ROI summary, not inferred from voxel totals."), "one target")
  fraction("target_unique", "percent_of_unique_focal_intersecting_streamlines", "unique_target_positive_focal_intersecting_streamlines", "unique_focal_intersecting_streamlines from the ROI summary", "one target")
  add("target_unique", "checksum", "Optional user-supplied target-mask checksum.", "checksum", "one target")
  counts("target_associations", c(focal_voxel_streamline_target_associations = "Sum of this target's per-voxel counts.", focal_voxels_with_at_least_one_association = "Focal voxels with a positive count for this target."), "one target")
  fraction("target_associations", "percent_of_all_focal_voxel_streamline_associations", "focal_voxel_streamline_target_associations", "sum of per-voxel tracks", "one target")
  for (nm in c("target", "axis")) add("spatial_gradients", nm, if (nm == "target") "Target label." else "Image array axis x, y or z, not a scanner-space anatomical direction.", population = "one target/axis")
  add("spatial_gradients", "slope", "Coefficient of the selected array axis in target_count/tracks ~ x + y + z, fitted over finite fractions at positive-assignment focal voxels.", "fraction per voxel index", "one target/axis", missing_zero = "NA when the coefficient is not estimable; unqualified targets have no rows.")
  add("spatial_gradients", "p_value", "Nominal linear-model coefficient p-value. Within-subject voxels are not independent biological replicates; this is not a group-level inference.", "probability", "one target/axis", missing_zero = "NA for a non-estimable coefficient.")
  add("spatial_gradients", "q_value", "Benjamini-Hochberg adjustment of p_value across returned target-by-axis rows.", "adjusted probability", "all returned target/axis rows", missing_zero = "NA when p_value is unavailable.")
  counts("spatial_gradients", c(connected_voxels = "Positive-assignment focal voxels with a finite target fraction included in the model.", target_connected_voxels = "Included voxels with a positive target count."), "one target model")

  # Optional additive statistical endpoint tables inherit the same endpoint
  # definitions and append labels, without changing their counting population.
  copy("endpoint_summary", "stat_endpoints")
  copy("stat_voxels", "stat_endpoints", c("subject", "group", "roi"))
  copy("penetration_summary", "stat_penetration")
  copy("stat_voxels", "stat_penetration", c("subject", "group", "roi"))
  copy("endpoint_summary", "target_endpoints")
  copy("target_unique", "target_endpoints", c("target_family", "checksum"))

  out <- do.call(rbind, rows)
  aliases <- c(
    "rsc/focal_voxel_streamline_target_associations" = "target_associations", "rsc/tracks" = "target_associations", "rsc/total_origin_voxels" = "total_focal_voxels",
    "roi_summary/unique_nonlocal_target_positive_streamlines" = "unique_leave_with_target_streamlines",
    "roi_summary/unique_nonlocal_target_negative_streamlines" = "unique_leave_no_target_streamlines",
    "roi_summary/unique_any_target_positive_streamlines" = "unique_focal_intersecting_with_target_streamlines",
    "penetration_summary/at_least_ten_voxels_streamlines" = "ten_or_more_voxels_streamlines",
    "penetration_summary/at_least_ten_voxels_percent_of_target_intersections" = "ten_or_more_voxels_percent_of_target_intersections"
  )
  keys <- paste(out$family, out$column, sep = "/")
  out$alias_of[match(names(aliases), keys)] <- unname(aliases)
  out$method[out$family %in% c("endpoint_summary", "endpoint_penetration", "penetration_summary", "stat_endpoints", "stat_penetration", "target_endpoints")] <- endpoint_method
  out$definition[out$family %in% c("endpoint_summary", "stat_endpoints", "target_endpoints") & out$column == "target_intersection_mode"] <- "Always full_streamline in this additive audit, even when the run's primary assignment uses endpoints or local search."
  out$definition[out$family == "run_index" & out$column == "run_status"] <- "Legacy aggregate state: completed only if the analysis and all requested stages completed without error; failed otherwise. Inspect analysis_status and per-stage fields to identify preserved results."
  out$missing_zero[out$family == "endpoint_context" & out$column %in% c("endpoints_nonfinite", "endpoints_one_voxel_outside", "endpoints_more_than_one_voxel_outside")] <- "0 means none; NA in anatomical-mask rows because this is a reference-grid diagnostic."
  rownames(out) <- NULL
  out
}

ftq_output_field_definition <- function(family, columns, target_labels = character(),
                                        coordinate_columns = character(),
                                        derived_column = character()) {
  dictionary <- ftq_output_field_dictionary()
  dictionary <- dictionary[dictionary$family == family, , drop = FALSE]
  if (!nrow(dictionary)) stop("Unknown output family: ", family, call. = FALSE)
  keys <- as.character(columns)
  if (family %in% c("per_voxel", "all_voxels", "top_voxels")) keys[keys %in% target_labels] <- "{target}"
  if (family %in% c("all_voxels", "top_voxels")) {
    coordinate_keys <- columns %in% coordinate_columns
    if (family == "top_voxels" && any(coordinate_keys)) {
      if (length(coordinate_columns) != 3L) stop("Exactly three ordered coordinate columns are required.", call. = FALSE)
      keys[coordinate_keys] <- c("x", "y", "z")[match(columns[coordinate_keys], coordinate_columns)]
    } else keys[coordinate_keys] <- "{coordinate}"
  }
  if (family == "top_voxels") keys[columns %in% derived_column] <- "{derived_leave_with_target}"
  if (family == "roi_summary") {
    for (prefix in c("unique_target_streamlines:", "unique_full_path_target_streamlines:")) {
      keys[columns %in% paste0(prefix, target_labels)] <- paste0(prefix, "{target}")
    }
  }
  at <- match(keys, dictionary$column)
  if (anyNA(at)) stop("Undocumented ", family, " column(s): ", paste(columns[is.na(at)], collapse = ", "), call. = FALSE)
  out <- dictionary[at, , drop = FALSE]
  out$column <- as.character(columns)
  rownames(out) <- NULL
  out
}

ftq_output_family_files <- function() {
  c(
    per_voxel = "*_fibrtraq_per-voxel-hits.csv",
    all_voxels = "tables/*_all-voxels.csv",
    all_voxels_columns = "tables/*_all-voxels-columns.csv",
    local_fraction = "tables/*_local-fraction.csv",
    rsc = "*_fibrtraq_rsc-summary.csv",
    fate = "*_fibrtraq_fate-summary.csv (voxel-weighted)",
    roi_summary = "*_fibrtraq_roi-summary.csv",
    connection_summary = "tables/*_connection-summary.csv",
    report_fate = "tables/*_report-fate-summary.csv",
    top_voxels = "tables/*_ranked_*_voxels.csv (all focal voxels; historical top-N previews remain readable)",
    endpoint_summary = "*_fibrtraq_endpoint-summary.csv",
    endpoint_penetration = "*_fibrtraq_endpoint-penetration.csv",
    penetration_summary = "tables/*_endpoint-penetration-summary.csv",
    endpoint_context = "*_fibrtraq_endpoint-context.csv (Endpoint coverage)",
    local_search_assignments = "*_fibrtraq_local-search-assignments.csv",
    local_search_summary = "tables/*_local-search-summary.csv",
    local_search_endpoint_status = "tables/*_local-search-endpoint-status.csv",
    exploratory_per_voxel = "ftq_exploratory_summary()$data; explicit research export: *_fibrtraq_exploratory-per-voxel.csv",
    exploratory_target_totals = "ftq_exploratory_summary()$target_totals; explicit research export: *_fibrtraq_exploratory-target-totals.csv",
    exploratory_lut = "*_fibrtraq_exploratory-primary-target-lut.csv",
    exploratory_provenance = "ftq_exploratory_summary()$provenance; explicit research export: *_fibrtraq_exploratory-provenance.csv",
    run_index = "fibrtraq-run-index.csv",
    prepared_inputs = "fibrtraq-prepared-inputs.csv",
    prepared_volumes = "prepared/*_prepared-reference-volumes.csv or *_prepared-background-volumes.csv",
    union_inputs = "prepared/*_regional-union-inputs.csv",
    union_overlaps = "prepared/*_regional-union-overlaps.csv",
    tractogram_context = "*_metadata.json: files.context.summary; explicit research export: *_fibrtraq_tractogram-context-summary.csv",
    mask_grid = "ftq_mask_grid_report(..., output_csv = ...) destination",
    validation_suite = "ftq_run_validation_suite(..., output_csv = ...) destination",
    stat_roi = "ftq_stat_subject_summary()$roi; saved by the installed study recipe",
    stat_targets = "ftq_stat_subject_summary()$targets; saved by the installed study recipe",
    stat_voxels = "ftq_stat_subject_summary()$voxels; saved by the installed study recipe",
    target_unique = "ftq_target_summary_tables()$unique_streamlines",
    target_associations = "ftq_target_summary_tables()$voxel_associations",
    spatial_gradients = "ftq_stat_spatial_gradients() return table; saved by the installed study recipe",
    stat_endpoints = "ftq_stat_subject_summary()$endpoints (optional)",
    stat_penetration = "ftq_stat_subject_summary()$penetration (optional)",
    target_endpoints = "ftq_target_summary_tables()$endpoint_streamlines (optional)"
  )
}

ftq_output_guide_sections <- function() {
  list(
    "Start with Results" = c(
      "Open the study or single-run report index and choose Results. Read the subject, focal region, assignment method, computational status and biological QC separately. The RSC summary, core maps and complete voxel table are the main results; advanced tables remain available for audit.",
      "The complete original *_fibrtraq_per-voxel-hits.csv contains every focal voxel. The readable *_all-voxels.csv adds explicit array coordinates; use *_all-voxels-columns.csv to identify coordinate and target roles safely. Ranked CSVs contain every focal voxel, including zero-count voxels: descending selected metric, then increasing x, y, z for ties. Historical top-50 files contain only a preview.",
      "Use the original per-voxel-hits CSV as the general analysis input. The all-voxel file is a readable export: a plain read.csv() retains coordinate columns, which must not be passed to target-inference functions as if they were target counts. FibrTraQ file readers recognize the named all-voxel CSV together with its companion column-role record and restore the canonical fields.",
      "Copy the complete result folder and all subfolders, then extract it before opening HTML. HTML uses portable relative links; recorded source paths in CSV/JSON preserve provenance from the generating workstation. JSON is machine-readable text, not PDF; CSV is a table, not an Excel workbook."
    ),
    "Report locations" = c(
      "The main Results report is under reports/<run>/; its images and embedded viewer are in report-figures/. Complete ranked CSVs are in tables/ beside the report. The supporting Target report and its five PNG figures are in the focal-region figures/ folder. The primary-target lookup CSV is in the focal-region folder. Standard runs do not create exploratory-report/ or context/ folders.",
      "The Target report retains its calculated summaries in memory and displays the method and denominators. For further research, ftq_exploratory_summary() can explicitly export the detailed derived tables and maps to a separate directory. Existing research exports and historical result folders remain readable."
    ),
    "Counting populations and missing values" = c(
      "A unique streamline count counts a tractogram record once within its stated population. A focal-voxel/streamline association counts its allocation to a particular focal voxel; summing across voxels can count a streamline repeatedly. Target-association totals can also repeat across targets. Never infer biological sample size from voxels, streamlines or workflow runs.",
      "RSC (%) = 100 * focal voxels with a positive count for this target / ALL focal voxels, including zero-assignment voxels. It is not a percentage of streamlines. A focal voxel can contribute to several targets, even when each streamline is assigned to only one target, so RSC percentages do not need to sum to 100%. Local means every selected path voxel lies in the configured local ROI, which may differ from the focal ROI. Local and target-positive connection rows can overlap. The three fate rows are exclusive because local takes precedence, and the two nonlocal rows split by target assignment.",
      "All focal voxels are retained. Under all_focal_voxels, a streamline is allocated to each distinct occupied focal voxel. Under first_intersection, only its first selected focal voxel receives the count: a zero-assignment voxel can still be traversed. Do not equate zero-assignment coverage with biological absence or endpoint coverage.",
      "NA means unavailable, inapplicable or unevaluable as specified for each field, not zero. Several legacy fraction fields store 0 for zero denominators: the dictionary identifies them. Exclude unevaluable voxels when interpreting distributions. The additive local-fraction CSV uses NA for zero tracks; its NIfTI stores 0 with a separate validity mask."
    ),
    "Method and endpoint audit" = c(
      "full_streamline assigns targets from selected path voxels; endpoint_voxels requires an unambiguous focal/target endpoint pair; local_search requires an unambiguous nearest-label pair within the inclusive physical radius. The complete target set affects endpoint/local-search assignment through overlap, ties and nearer labels. Preserve that set and the recorded method when rebuilding reports.",
      "The additive endpoint-summary and penetration tables always audit exact full-path target intersection. Exact endpoint membership and dilated1 (Chebyshev radius one voxel, including 26 neighbours) are separate from physical-radius local search. Broad paired_endpoint counts may include shared-mask ambiguity; the unambiguous paired_focal_target_endpoints class excludes it. Use the numerator named in each percentage definition.",
      "Endpoint coverage counts all valid mapped streamlines before focal filtering. 100% inside the image grid means endpoints lie inside image bounds; it does not establish anatomy, continuous path coverage, local classification or target assignment. With trusted masks, anatomical endpoint containment remains QC evidence rather than proof of biological validity.",
      "Penetration is the number of unique exact target voxels occupied by the selected full path. It is not physical depth. Histogram bins are sparse; their counts sum to full-path target-intersection streamlines. Under sample_points, sparse stored points can miss intervening voxels. The five exclusive endpoint classes likewise partition that same target-intersection population."
    ),
    "Maps, coordinates and image lookup" = c(
      "fNode_Coords and readable x/y/z values are one-based image-array indices, not scanner millimetres. Image affines and declared spatial units supply the physical transformation. Existing target names, including x, y or z, are preserved; the companion column-role table identifies renamed coordinate fields.",
      "Core maps in maps/ retain the focal reference grid and affine: tracks (assigned streamline count), tracks-binary (1 if tracks > 0), local, leave-no-target, leave-with-target, and target-{label}. Outside the focal mask these maps are zero. The original focal mask is necessary to distinguish outside-mask zeros from zero-assignment focal voxels. Display scaling affects previews only.",
      "The additive local-fraction map is local / tracks on positive-assignment focal voxels, with stored zero elsewhere. The accompanying local-fraction-valid map is 1 only where tracks > 0, including zero-fraction but evaluable voxels. Read them together. These maps use the configured local ROI, not an assumed focal-only region.",
      "Exploratory maps are available only through the explicit research export ftq_exploratory_summary(..., output_dir = ..., template_image = ..., write_maps = TRUE). Standard reports do not write these maps or their supporting CSVs. Research maps retain the corresponding exploratory-per-voxel quantities: primary_target_index, primary_target_hits, primary_target_fraction, target_diversity, target_entropy_bits, target_entropy_fraction, total_target_hits, target_hits_per_assigned_track, leave_no_target_fraction and leave_with_target_fraction. They are supporting views. Category index 0 means no target; positive indices use the primary-target LUT. Ties choose the first target column. A stored zero entropy can mean no data or a concentrated distribution, so inspect target/assignment coverage.",
      "Tractogram density context is whole-tractogram streamline occupancy per reference voxel (at most once per streamline per voxel). It is a context image, not MRI, a brain mask or anatomical tissue volume. Core maps and context maps have different counting populations. The density image is saved with the core maps; its summary and settings are retained in the run metadata. PNG figures are previews; their portable download links provide full-size images, while NIfTI files preserve the full voxel values.",
      "The 3-D viewer is a quota-limited, scan-order subset selected by target assignment. It is neither random nor prevalence weighted. Its paths and endpoint markers support spatial inspection; they must not be used to estimate target prevalence or whole-tractogram completeness."
    ),
    "Methods, validation and supporting records" = c(
      "*_metadata.json preserves generation time, loaded package identity, normalized manifest, input paths/checksums, compact preflight, counting summaries, endpoint method, local-search geometry and output paths. *_validation.json preserves computational check results, reconciliation and expected-baseline checks where requested, plus separate biological QC warnings/status. Readable Methods and Validation sections expose these records without replacing their JSON source.",
      "Prepared input, union and selected-volume CSVs document what was used or derived. The prepared/ folder is created only when an image or mask derivative is needed. Optional unavailable inputs are not analysis failures. Run indexes identify every subject/focal analysis and its outputs; record and inspect analysis versus presentation status separately. Local output paths and subject labels may contain information that needs review before sharing publicly.",
      "The bundled synthetic example compares seven scientific tables and 18 independently declared claims exactly. Four HTML artifacts are checked for existence and nonzero size. These checks do not compare rendered HTML exactly, exercise browser interaction or establish biological validity.",
      "The families below describe standard file-workflow outputs, explicitly requested research exports, historical fields and the public grid/validation diagnostic CSVs. Target labels and added coordinates are explicit parameterized roles, not inferred numeric columns. Optional columns/whole tables appear only when their source feature is enabled. All legacy aliases remain documented so saved analyses remain readable."
    )
  )
}

ftq_output_guide_markdown <- function() {
  lines <- c("# FibrTraQ output guide", "", "Definitions for the installed package's public file workflows. Keep this guide with the complete results folder.", "")
  sections <- ftq_output_guide_sections()
  for (nm in names(sections)) lines <- c(lines, paste0("## ", nm), "", as.vector(rbind(sections[[nm]], "")))
  dictionary <- ftq_output_field_dictionary()
  files <- ftq_output_family_files()
  for (family in unique(dictionary$family)) {
    part <- dictionary[dictionary$family == family, , drop = FALSE]
    lines <- c(lines, paste0("## ", family), "", paste0("Files or returned table: `", files[[family]], "`."), "")
    for (i in seq_len(nrow(part))) {
      row <- part[i, ]
      lines <- c(lines, paste0("### `", row$column, "`"), "", row$definition, "",
                 paste0("- Units: ", row$units, "."), paste0("- Population: ", row$population, "."),
                 paste0("- Denominator: ", row$denominator, "."), paste0("- Missing/zero: ", row$missing_zero),
                 paste0("- Method: ", row$method, "."), paste0("- Alias of: ", row$alias_of, "."), "")
    }
  }
  lines[-length(lines)]
}

ftq_output_guide_html <- function(markdown_file = "OUTPUTS_README.md") {
  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  dictionary <- ftq_output_field_dictionary()
  files <- ftq_output_family_files()
  families <- unique(dictionary$family)
  html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1"><title>FibrTraQ output guide</title>',
    '<style>body{font:17px/1.55 system-ui,sans-serif;color:#182735;background:#f6f8fa;margin:0}main{max-width:1000px;margin:auto;padding:2rem}h1,h2{line-height:1.2}a{color:#075a9d}details{background:white;border:1px solid #cbd5df;border-radius:.4rem;margin:1rem 0;padding:1rem}summary{font-weight:650;cursor:pointer;overflow-wrap:anywhere}dt{margin-top:1.5rem;font-weight:650;overflow-wrap:anywhere}dd{margin:.4rem 0 1.5rem}dd p{margin:.35rem 0}code{font-size:.9em}nav ul{columns:2;padding-left:1.3rem}li{break-inside:avoid;overflow-wrap:anywhere}.meta{font-size:.9rem;color:#3c4e5e}.note{border-left:4px solid #32759b;padding-left:1rem}@media(max-width:600px){main{padding:1rem}nav ul{columns:1}}@media print{details{break-inside:auto}summary{font-size:1.2rem}}</style></head><body><main>',
    '<h1>FibrTraQ output guide</h1><p>Definitions for the installed package\u2019s public file workflows. Keep this guide with the complete results folder.</p>',
    paste0('<p><a href="', esc(markdown_file), '">Download Markdown definitions</a></p>'))
  sections <- ftq_output_guide_sections()
  for (nm in names(sections)) html <- c(html, paste0("<h2>", esc(nm), "</h2>"), paste0("<p>", esc(sections[[nm]]), "</p>"))
  html <- c(html, '<h2>Find a table family</h2><nav aria-label="Output families"><ul>',
            paste0('<li><a href="#', families, '">', gsub("_", " ", families), '</a></li>'), '</ul></nav>')
  for (family in families) {
    part <- dictionary[dictionary$family == family, , drop = FALSE]
    html <- c(html, paste0('<details id="', family, '"><summary>', gsub("_", " ", family), ' \u2014 ', nrow(part), ' field definitions</summary><p><strong>Files or returned table:</strong> <code>', esc(files[[family]]), '</code></p><dl>'))
    for (i in seq_len(nrow(part))) {
      row <- part[i, ]
      html <- c(html, paste0('<dt><code>', esc(row$column), '</code></dt><dd><p>', esc(row$definition), '</p>'),
                paste0('<p class="meta"><strong>Units:</strong> ', esc(row$units), '; <strong>Population:</strong> ', esc(row$population),
                       '; <strong>Denominator:</strong> ', esc(row$denominator), '.</p>'),
                paste0('<p class="meta"><strong>Missing/zero:</strong> ', esc(row$missing_zero), '</p>'),
                paste0('<p class="meta"><strong>Method:</strong> ', esc(row$method), '; <strong>Alias of:</strong> ', esc(row$alias_of), '.</p></dd>'))
    }
    html <- c(html, '</dl></details>')
  }
  c(html, '</main></body></html>')
}

ftq_write_output_guide <- function(output_dir) {
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one nonempty directory path.", call. = FALSE)
  }
  output_dir <- path.expand(output_dir)
  contents <- list(OUTPUTS_README.html = ftq_output_guide_html(), OUTPUTS_README.md = ftq_output_guide_markdown())
  paths <- file.path(output_dir, names(contents))
  bytes <- function(lines) charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n")))
  expected <- lapply(contents, bytes)
  differs <- function(path, value) {
    if (!file.exists(path)) return(FALSE)
    size <- file.info(path)$size
    dir.exists(path) || !is.finite(size) || size != length(value) ||
      !identical(readBin(path, "raw", n = size), value)
  }
  # A new package version can render into an existing results directory.
  # Preserve customized/older guides and select a stable content-specific pair.
  if (any(vapply(seq_along(paths), function(i) differs(paths[[i]], expected[[i]]), logical(1)))) {
    fingerprint_file <- tempfile("fibrtraq-guide-")
    on.exit(unlink(fingerprint_file), add = TRUE)
    writeBin(c(expected[[1L]], expected[[2L]]), fingerprint_file)
    suffix <- substr(unname(tools::md5sum(fingerprint_file)), 1L, 12L)
    names(contents) <- paste0("OUTPUTS_README-", suffix, c(".html", ".md"))
    contents[[1L]] <- ftq_output_guide_html(names(contents)[[2L]])
    paths <- file.path(output_dir, names(contents))
    expected <- lapply(contents, bytes)
  }
  # Check both destinations before writing either, including the uncommon case
  # of a user modifying a previously versioned guide in place.
  for (i in seq_along(paths)) if (differs(paths[[i]], expected[[i]])) {
    stop("Content-specific output guide already exists with different content: ", paths[[i]],
         ". Keep it and choose a new output directory.", call. = FALSE)
  }
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output guide directory: ", output_dir, call. = FALSE)
  }
  for (i in seq_along(paths)) if (!file.exists(paths[[i]])) writeBin(expected[[i]], paths[[i]])
  paths[[1L]]
}
