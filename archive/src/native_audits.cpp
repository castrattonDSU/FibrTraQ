#include "../../src/voxel_path.h"
#include <unordered_set>

using namespace Rcpp;

namespace {
void ensure_hist_size(std::vector<double> &hist, double value, double bin_width) {
  if (!std::isfinite(value) || value < 0.0) return;
  const std::size_t bin = static_cast<std::size_t>(std::floor(value / bin_width));
  if (hist.size() <= bin) hist.resize(bin + 1, 0.0);
  hist[bin] += 1.0;
}
} // namespace

// [[Rcpp::export]]
List ftq_spacing_chunk_cpp(List tracks,
                           NumericMatrix raw_to_voxel,
                           NumericVector voxel_size,
                           double bin_width_mm = 0.01,
                           double bin_width_vox = 0.01) {
  if (raw_to_voxel.nrow() < 3 || raw_to_voxel.ncol() < 4) {
    stop("raw_to_voxel must have at least 3 rows and 4 columns.");
  }
  if (voxel_size.size() < 3) {
    stop("voxel_size must have length at least 3.");
  }
  if (bin_width_mm <= 0.0 || bin_width_vox <= 0.0) {
    stop("bin widths must be positive.");
  }

  double segments = 0.0;
  double axis_gt_one = 0.0;
  double euclidean_gt_one_edge = 0.0;
  double euclidean_gt_one_diagonal = 0.0;
  double min_mm = R_PosInf;
  double max_mm = R_NegInf;
  double min_vox = R_PosInf;
  double max_vox = R_NegInf;
  double min_dx = R_PosInf;
  double max_dx = R_NegInf;
  double min_dy = R_PosInf;
  double max_dy = R_NegInf;
  double min_dz = R_PosInf;
  double max_dz = R_NegInf;
  std::vector<double> hist_mm;
  std::vector<double> hist_vox;
  std::vector<double> hist_dx;
  std::vector<double> hist_dy;
  std::vector<double> hist_dz;

  for (R_xlen_t t = 0; t < tracks.size(); ++t) {
    NumericMatrix coords = track_coords(tracks[t]);
    if (coords.ncol() < 3 || coords.nrow() < 2) continue;

    double prev_x;
    double prev_y;
    double prev_z;
    raw_to_voxel_point(coords(0, 0), coords(0, 1), coords(0, 2),
                       raw_to_voxel, prev_x, prev_y, prev_z);
    for (int i = 1; i < coords.nrow(); ++i) {
      double vx;
      double vy;
      double vz;
      raw_to_voxel_point(coords(i, 0), coords(i, 1), coords(i, 2),
                         raw_to_voxel, vx, vy, vz);
      if (std::isfinite(prev_x) && std::isfinite(prev_y) && std::isfinite(prev_z) &&
          std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz)) {
        const double dx = std::fabs(vx - prev_x);
        const double dy = std::fabs(vy - prev_y);
        const double dz = std::fabs(vz - prev_z);
        const double dmm_x = dx * voxel_size[0];
        const double dmm_y = dy * voxel_size[1];
        const double dmm_z = dz * voxel_size[2];
        const double euclidean_mm = std::sqrt(dmm_x * dmm_x + dmm_y * dmm_y + dmm_z * dmm_z);
        const double euclidean_vox = std::sqrt(dx * dx + dy * dy + dz * dz);

        segments += 1.0;
        if (euclidean_mm < min_mm) min_mm = euclidean_mm;
        if (euclidean_mm > max_mm) max_mm = euclidean_mm;
        if (euclidean_vox < min_vox) min_vox = euclidean_vox;
        if (euclidean_vox > max_vox) max_vox = euclidean_vox;
        if (dx < min_dx) min_dx = dx;
        if (dx > max_dx) max_dx = dx;
        if (dy < min_dy) min_dy = dy;
        if (dy > max_dy) max_dy = dy;
        if (dz < min_dz) min_dz = dz;
        if (dz > max_dz) max_dz = dz;
        if (dx > 1.0 || dy > 1.0 || dz > 1.0) axis_gt_one += 1.0;
        if (euclidean_vox > 1.0) euclidean_gt_one_edge += 1.0;
        if (euclidean_vox > std::sqrt(3.0)) euclidean_gt_one_diagonal += 1.0;
        ensure_hist_size(hist_mm, euclidean_mm, bin_width_mm);
        ensure_hist_size(hist_vox, euclidean_vox, bin_width_vox);
        ensure_hist_size(hist_dx, dx, bin_width_vox);
        ensure_hist_size(hist_dy, dy, bin_width_vox);
        ensure_hist_size(hist_dz, dz, bin_width_vox);
      }
      prev_x = vx;
      prev_y = vy;
      prev_z = vz;
    }
  }

  return List::create(
    _["streamlines"] = static_cast<double>(tracks.size()),
    _["segments"] = segments,
    _["min_mm"] = min_mm,
    _["max_mm"] = max_mm,
    _["min_vox"] = min_vox,
    _["max_vox"] = max_vox,
    _["min_dx"] = min_dx,
    _["max_dx"] = max_dx,
    _["min_dy"] = min_dy,
    _["max_dy"] = max_dy,
    _["min_dz"] = min_dz,
    _["max_dz"] = max_dz,
    _["axis_gt_one"] = axis_gt_one,
    _["euclidean_gt_one_edge"] = euclidean_gt_one_edge,
    _["euclidean_gt_one_diagonal"] = euclidean_gt_one_diagonal,
    _["hist_mm"] = NumericVector(hist_mm.begin(), hist_mm.end()),
    _["hist_vox"] = NumericVector(hist_vox.begin(), hist_vox.end()),
    _["hist_dx"] = NumericVector(hist_dx.begin(), hist_dx.end()),
    _["hist_dy"] = NumericVector(hist_dy.begin(), hist_dy.end()),
    _["hist_dz"] = NumericVector(hist_dz.begin(), hist_dz.end())
  );
}

// [[Rcpp::export]]
List ftq_endpoint_audit_chunk_cpp(List tracks,
                                  NumericMatrix raw_to_voxel,
                                  std::string method,
                                  LogicalVector axis_flip,
                                  IntegerVector vol_dim,
                                  LogicalMatrix region_lut,
                                  LogicalMatrix region_lut_dilated,
                                  LogicalMatrix context_lut,
                                  std::string streamline_voxelization_mode,
                                  int penetration_hist_max = 500) {
  if (raw_to_voxel.nrow() < 3 || raw_to_voxel.ncol() < 4) {
    stop("raw_to_voxel must have at least 3 rows and 4 columns.");
  }
  if (axis_flip.size() != 3 || vol_dim.size() < 3) {
    stop("axis_flip and vol_dim must have length 3.");
  }
  const int n_voxels = vol_dim[0] * vol_dim[1] * vol_dim[2];
  if (region_lut.nrow() != n_voxels ||
      region_lut_dilated.nrow() != n_voxels ||
      context_lut.nrow() != n_voxels) {
    stop("All lookup matrices must have one row per reference voxel.");
  }
  if (region_lut.ncol() != region_lut_dilated.ncol()) {
    stop("Exact and dilated region lookup matrices must have matching columns.");
  }
  if (penetration_hist_max < 1) {
    stop("penetration_hist_max must be positive.");
  }

  const std::vector<int> codes = method_codes(CharacterVector::create(method));
  const int method_code = codes[0];
  const bool interpolate_segments =
    streamline_voxelization_mode == "segment_interpolation";
  if (streamline_voxelization_mode != "sample_points" && !interpolate_segments) {
    stop("Unsupported streamline_voxelization_mode in endpoint audit: " +
         streamline_voxelization_mode);
  }

  const int n_regions = region_lut.ncol();
  const int n_contexts = context_lut.ncol();
  IntegerVector region_intersection(n_regions);
  IntegerVector endpoint_region_exact(n_regions);
  IntegerVector endpoint_region_dilated(n_regions);
  IntegerMatrix intersection_pair(n_regions, n_regions);
  IntegerMatrix intersection_pair_target_endpoint_exact(n_regions, n_regions);
  IntegerMatrix intersection_pair_target_endpoint_dilated(n_regions, n_regions);
  IntegerMatrix intersection_pair_focal_endpoint_exact(n_regions, n_regions);
  IntegerMatrix intersection_pair_focal_endpoint_dilated(n_regions, n_regions);
  IntegerMatrix intersection_pair_any_relevant_endpoint_exact(n_regions, n_regions);
  IntegerMatrix intersection_pair_any_relevant_endpoint_dilated(n_regions, n_regions);
  IntegerMatrix focal_endpoint_target_intersection_exact(n_regions, n_regions);
  IntegerMatrix focal_endpoint_target_intersection_dilated(n_regions, n_regions);
  IntegerMatrix endpoint_pair_exact(n_regions, n_regions);
  IntegerMatrix endpoint_pair_dilated(n_regions, n_regions);
  IntegerMatrix intersection_pair_endpoint_pair_dilated(n_regions, n_regions);
  IntegerMatrix endpoint_class_exact(n_regions * n_regions, 5);
  IntegerMatrix endpoint_class_dilated(n_regions * n_regions, 5);
  IntegerMatrix penetration_hist(
    n_regions * n_regions,
    penetration_hist_max + 1
  );
  IntegerVector endpoint_context(n_contexts);
  IntegerVector streamline_any_endpoint_context(n_contexts);
  IntegerVector streamline_both_endpoints_context(n_contexts);
  IntegerVector endpoint_out_of_bounds_low(3);
  IntegerVector endpoint_out_of_bounds_high(3);
  IntegerMatrix endpoint_out_of_bounds_low_distance(3, 2);
  IntegerMatrix endpoint_out_of_bounds_high_distance(3, 2);

  int processed = 0;
  int valid_streamlines = 0;
  int endpoint_count = 0;
  int endpoint_out_of_bounds = 0;

  for (R_xlen_t t = 0; t < tracks.size(); ++t) {
    NumericMatrix coords = track_coords(tracks[t]);
    ++processed;
    if (coords.ncol() < 3 || coords.nrow() == 0) continue;
    ++valid_streamlines;

    std::vector<int> keys = voxelized_track_keys(
      coords, raw_to_voxel, method_code, axis_flip, vol_dim,
      interpolate_segments
    );
    std::vector<int> intersects(n_regions, 0);
    std::vector<std::unordered_set<int> > region_voxels(n_regions);
    std::unordered_set<int> unique_track_voxels;
    unique_track_voxels.reserve(keys.size());
    for (std::size_t k = 0; k < keys.size(); ++k) {
      const int key = keys[k];
      if (key < 0 || !unique_track_voxels.insert(key).second) continue;
      for (int region = 0; region < n_regions; ++region) {
        if (region_lut(key, region) == TRUE) {
          intersects[region] = 1;
          region_voxels[region].insert(key);
        }
      }
    }

    int start_key = -1;
    int end_key = -1;
    discretized_key(
      coords(0, 0), coords(0, 1), coords(0, 2), raw_to_voxel,
      method_code, axis_flip, vol_dim, start_key
    );
    if (coords.nrow() == 1) {
      end_key = start_key;
    } else {
      discretized_key(
        coords(coords.nrow() - 1, 0),
        coords(coords.nrow() - 1, 1),
        coords(coords.nrow() - 1, 2),
        raw_to_voxel, method_code, axis_flip, vol_dim, end_key
      );
    }

    const int n_endpoints = coords.nrow() == 1 ? 1 : 2;
    endpoint_count += n_endpoints;
    if (start_key < 0) ++endpoint_out_of_bounds;
    if (n_endpoints == 2 && end_key < 0) ++endpoint_out_of_bounds;

    const auto record_endpoint_bounds = [&](int row) {
      double vx;
      double vy;
      double vz;
      raw_to_voxel_point(
        coords(row, 0), coords(row, 1), coords(row, 2),
        raw_to_voxel, vx, vy, vz
      );
      const double values[3] = {vx, vy, vz};
      for (int axis = 0; axis < 3; ++axis) {
        if (!std::isfinite(values[axis])) continue;
        double discrete;
        if (method_code == 0) {
          discrete = std::nearbyint(values[axis]);
        } else if (method_code == 1) {
          discrete = std::floor(values[axis]);
        } else {
          discrete = std::ceil(values[axis]);
        }
        int coordinate = static_cast<int>(discrete);
        if (axis_flip[axis] == TRUE) {
          coordinate = vol_dim[axis] - coordinate + 1;
        }
        if (coordinate < 1) {
          ++endpoint_out_of_bounds_low[axis];
          const int distance = 1 - coordinate;
          ++endpoint_out_of_bounds_low_distance(axis, distance == 1 ? 0 : 1);
        }
        if (coordinate > vol_dim[axis]) {
          ++endpoint_out_of_bounds_high[axis];
          const int distance = coordinate - vol_dim[axis];
          ++endpoint_out_of_bounds_high_distance(axis, distance == 1 ? 0 : 1);
        }
      }
    };
    record_endpoint_bounds(0);
    if (n_endpoints == 2) record_endpoint_bounds(coords.nrow() - 1);

    std::vector<int> start_exact(n_regions, 0);
    std::vector<int> end_exact(n_regions, 0);
    std::vector<int> start_dilated(n_regions, 0);
    std::vector<int> end_dilated(n_regions, 0);
    for (int region = 0; region < n_regions; ++region) {
      if (start_key >= 0) {
        start_exact[region] = region_lut(start_key, region) == TRUE;
        start_dilated[region] = region_lut_dilated(start_key, region) == TRUE;
      }
      if (end_key >= 0) {
        end_exact[region] = region_lut(end_key, region) == TRUE;
        end_dilated[region] = region_lut_dilated(end_key, region) == TRUE;
      }
    }

    std::vector<int> endpoint_any_exact(n_regions, 0);
    std::vector<int> endpoint_any_dilated(n_regions, 0);
    for (int region = 0; region < n_regions; ++region) {
      endpoint_any_exact[region] = start_exact[region] || end_exact[region];
      endpoint_any_dilated[region] = start_dilated[region] || end_dilated[region];
      region_intersection[region] += intersects[region];
      endpoint_region_exact[region] += endpoint_any_exact[region];
      endpoint_region_dilated[region] += endpoint_any_dilated[region];
    }

    for (int context = 0; context < n_contexts; ++context) {
      const bool start_inside =
        start_key >= 0 && context_lut(start_key, context) == TRUE;
      const bool end_inside =
        end_key >= 0 && context_lut(end_key, context) == TRUE;
      endpoint_context[context] += static_cast<int>(start_inside);
      if (n_endpoints == 2) {
        endpoint_context[context] += static_cast<int>(end_inside);
      }
      streamline_any_endpoint_context[context] +=
        static_cast<int>(start_inside || end_inside);
      streamline_both_endpoints_context[context] +=
        static_cast<int>(start_inside && (n_endpoints == 1 || end_inside));
    }

    for (int focal = 0; focal < n_regions; ++focal) {
      for (int target = 0; target < n_regions; ++target) {
        const bool pair_intersects = intersects[focal] && intersects[target];
        const bool shared_exact =
          (start_exact[focal] && start_exact[target]) ||
          (end_exact[focal] && end_exact[target]);
        const bool shared_dilated =
          (start_dilated[focal] && start_dilated[target]) ||
          (end_dilated[focal] && end_dilated[target]);
        const bool paired_exact = n_endpoints == 2 &&
          ((start_exact[focal] && end_exact[target]) ||
           (end_exact[focal] && start_exact[target]));
        const bool paired_dilated = n_endpoints == 2 &&
          ((start_dilated[focal] && end_dilated[target]) ||
           (end_dilated[focal] && start_dilated[target]));
        const int class_exact = shared_exact ? 0 :
          (paired_exact ? 1 :
           (endpoint_any_exact[target] ? 2 :
            (endpoint_any_exact[focal] ? 3 : 4)));
        const int class_dilated = shared_dilated ? 0 :
          (paired_dilated ? 1 :
           (endpoint_any_dilated[target] ? 2 :
            (endpoint_any_dilated[focal] ? 3 : 4)));

        intersection_pair(focal, target) += static_cast<int>(pair_intersects);
        if (pair_intersects) {
          endpoint_class_exact(focal * n_regions + target, class_exact) += 1;
          endpoint_class_dilated(focal * n_regions + target, class_dilated) += 1;
          intersection_pair_target_endpoint_exact(focal, target) +=
            endpoint_any_exact[target];
          intersection_pair_target_endpoint_dilated(focal, target) +=
            endpoint_any_dilated[target];
          intersection_pair_focal_endpoint_exact(focal, target) +=
            endpoint_any_exact[focal];
          intersection_pair_focal_endpoint_dilated(focal, target) +=
            endpoint_any_dilated[focal];
          intersection_pair_any_relevant_endpoint_exact(focal, target) +=
            static_cast<int>(
              endpoint_any_exact[focal] || endpoint_any_exact[target]
            );
          intersection_pair_any_relevant_endpoint_dilated(focal, target) +=
            static_cast<int>(
              endpoint_any_dilated[focal] || endpoint_any_dilated[target]
            );
          intersection_pair_endpoint_pair_dilated(focal, target) +=
            static_cast<int>(paired_dilated);
          const int unique_target_voxels =
            static_cast<int>(region_voxels[target].size());
          const int bin = std::min(unique_target_voxels, penetration_hist_max);
          penetration_hist(focal * n_regions + target, bin) += 1;
        }
        focal_endpoint_target_intersection_exact(focal, target) +=
          static_cast<int>(endpoint_any_exact[focal] && intersects[target]);
        focal_endpoint_target_intersection_dilated(focal, target) +=
          static_cast<int>(endpoint_any_dilated[focal] && intersects[target]);
        endpoint_pair_exact(focal, target) += static_cast<int>(paired_exact);
        endpoint_pair_dilated(focal, target) += static_cast<int>(paired_dilated);
      }
    }
  }

  return List::create(
    _["processed"] = processed,
    _["valid_streamlines"] = valid_streamlines,
    _["endpoint_count"] = endpoint_count,
    _["endpoint_out_of_bounds"] = endpoint_out_of_bounds,
    _["endpoint_out_of_bounds_low"] = endpoint_out_of_bounds_low,
    _["endpoint_out_of_bounds_high"] = endpoint_out_of_bounds_high,
    _["endpoint_out_of_bounds_low_distance"] =
      endpoint_out_of_bounds_low_distance,
    _["endpoint_out_of_bounds_high_distance"] =
      endpoint_out_of_bounds_high_distance,
    _["region_intersection"] = region_intersection,
    _["endpoint_region_exact"] = endpoint_region_exact,
    _["endpoint_region_dilated"] = endpoint_region_dilated,
    _["intersection_pair"] = intersection_pair,
    _["intersection_pair_target_endpoint_exact"] =
      intersection_pair_target_endpoint_exact,
    _["intersection_pair_target_endpoint_dilated"] =
      intersection_pair_target_endpoint_dilated,
    _["intersection_pair_focal_endpoint_exact"] =
      intersection_pair_focal_endpoint_exact,
    _["intersection_pair_focal_endpoint_dilated"] =
      intersection_pair_focal_endpoint_dilated,
    _["intersection_pair_any_relevant_endpoint_exact"] =
      intersection_pair_any_relevant_endpoint_exact,
    _["intersection_pair_any_relevant_endpoint_dilated"] =
      intersection_pair_any_relevant_endpoint_dilated,
    _["focal_endpoint_target_intersection_exact"] =
      focal_endpoint_target_intersection_exact,
    _["focal_endpoint_target_intersection_dilated"] =
      focal_endpoint_target_intersection_dilated,
    _["endpoint_pair_exact"] = endpoint_pair_exact,
    _["endpoint_pair_dilated"] = endpoint_pair_dilated,
    _["intersection_pair_endpoint_pair_dilated"] =
      intersection_pair_endpoint_pair_dilated,
    _["endpoint_class_exact"] = endpoint_class_exact,
    _["endpoint_class_dilated"] = endpoint_class_dilated,
    _["penetration_hist"] = penetration_hist,
    _["endpoint_context"] = endpoint_context,
    _["streamline_any_endpoint_context"] = streamline_any_endpoint_context,
    _["streamline_both_endpoints_context"] = streamline_both_endpoints_context
  );
}
