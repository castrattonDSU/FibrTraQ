#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector ftq_dilate_lut_cpp(LogicalVector lut,
                                 IntegerVector vol_dim) {
  if (vol_dim.size() != 3) stop("vol_dim must have length 3.");
  const int nx = vol_dim[0];
  const int ny = vol_dim[1];
  const int nz = vol_dim[2];
  if (nx <= 0 || ny <= 0 || nz <= 0) {
    stop("vol_dim values must be positive.");
  }
  const R_xlen_t expected =
    static_cast<R_xlen_t>(nx) * ny * nz;
  if (lut.size() != expected) {
    stop("Endpoint LUT dimensions are inconsistent.");
  }
  LogicalVector out = clone(lut);
  const int nxy = nx * ny;
  for (R_xlen_t key = 0; key < lut.size(); ++key) {
    if (lut[key] != TRUE) continue;
    const int index = static_cast<int>(key);
    const int z = index / nxy;
    const int remainder = index - z * nxy;
    const int y = remainder / nx;
    const int x = remainder - y * nx;
    const int x_low = std::max(0, x - 1);
    const int x_high = std::min(nx - 1, x + 1);
    const int y_low = std::max(0, y - 1);
    const int y_high = std::min(ny - 1, y + 1);
    const int z_low = std::max(0, z - 1);
    const int z_high = std::min(nz - 1, z + 1);
    for (int dz = z_low; dz <= z_high; ++dz) {
      const int z_offset = dz * nxy;
      for (int dy = y_low; dy <= y_high; ++dy) {
        const int yz_offset = z_offset + dy * nx;
        for (int dx = x_low; dx <= x_high; ++dx) {
          out[yz_offset + dx] = TRUE;
        }
      }
    }
  }
  return out;
}

namespace {

int encode_key(int x, int y, int z, const IntegerVector &vol_dim) {
  const int nx = vol_dim[0];
  const int ny = vol_dim[1];
  const int nz = vol_dim[2];
  if (x < 1 || x > nx || y < 1 || y > ny || z < 1 || z > nz) return -1;
  return (x - 1) + (y - 1) * nx + (z - 1) * nx * ny;
}

bool voxel_key(double vx,
               double vy,
               double vz,
               int method,
               const LogicalVector &axis_flip,
               const IntegerVector &vol_dim,
               int &key) {
  if (!std::isfinite(vx) || !std::isfinite(vy) || !std::isfinite(vz)) {
    key = -1;
    return false;
  }

  double dx;
  double dy;
  double dz;
  if (method == 0) {
    // R's round() uses IEC 60559 round-to-even for .5 ties. nearbyint() uses
    // the process rounding mode, which R leaves at the default round-to-even.
    dx = std::nearbyint(vx);
    dy = std::nearbyint(vy);
    dz = std::nearbyint(vz);
  } else if (method == 1) {
    dx = std::floor(vx);
    dy = std::floor(vy);
    dz = std::floor(vz);
  } else {
    dx = std::ceil(vx);
    dy = std::ceil(vy);
    dz = std::ceil(vz);
  }

  int x = static_cast<int>(dx);
  int y = static_cast<int>(dy);
  int z = static_cast<int>(dz);
  if (axis_flip[0] == TRUE) x = vol_dim[0] - x + 1;
  if (axis_flip[1] == TRUE) y = vol_dim[1] - y + 1;
  if (axis_flip[2] == TRUE) z = vol_dim[2] - z + 1;

  key = encode_key(x, y, z, vol_dim);
  return key >= 0;
}

bool discretized_key(double raw_x,
                     double raw_y,
                     double raw_z,
                     const NumericMatrix &raw_to_voxel,
                     int method,
                     const LogicalVector &axis_flip,
                     const IntegerVector &vol_dim,
                     int &key) {
  const double vx = raw_to_voxel(0, 0) * raw_x +
    raw_to_voxel(0, 1) * raw_y +
    raw_to_voxel(0, 2) * raw_z +
    raw_to_voxel(0, 3);
  const double vy = raw_to_voxel(1, 0) * raw_x +
    raw_to_voxel(1, 1) * raw_y +
    raw_to_voxel(1, 2) * raw_z +
    raw_to_voxel(1, 3);
  const double vz = raw_to_voxel(2, 0) * raw_x +
    raw_to_voxel(2, 1) * raw_y +
    raw_to_voxel(2, 2) * raw_z +
    raw_to_voxel(2, 3);
  return voxel_key(vx, vy, vz, method, axis_flip, vol_dim, key);
}

NumericMatrix track_coords(SEXP track) {
  if (Rf_isMatrix(track)) {
    return as<NumericMatrix>(track);
  }
  if (TYPEOF(track) == VECSXP) {
    List item(track);
    if (item.containsElementNamed("coords")) {
      return as<NumericMatrix>(item["coords"]);
    }
  }
  stop("Each track must be a matrix or a list with a coords matrix.");
}

std::vector<int> method_codes(CharacterVector methods) {
  std::vector<int> out;
  out.reserve(methods.size());
  for (R_xlen_t i = 0; i < methods.size(); ++i) {
    const std::string method = as<std::string>(methods[i]);
    if (method == "round") {
      out.push_back(0);
    } else if (method == "floor") {
      out.push_back(1);
    } else if (method == "ceil") {
      out.push_back(2);
    } else {
      stop("Unsupported discretization method in C++ backend: " + method);
    }
  }
  return out;
}

void ensure_hist_size(std::vector<double> &hist, double value, double bin_width) {
  if (!std::isfinite(value) || value < 0.0) return;
  const std::size_t bin = static_cast<std::size_t>(std::floor(value / bin_width));
  if (hist.size() <= bin) hist.resize(bin + 1, 0.0);
  hist[bin] += 1.0;
}

void raw_to_voxel_point(double raw_x,
                        double raw_y,
                        double raw_z,
                        const NumericMatrix &raw_to_voxel,
                        double &vx,
                        double &vy,
                        double &vz) {
  vx = raw_to_voxel(0, 0) * raw_x +
    raw_to_voxel(0, 1) * raw_y +
    raw_to_voxel(0, 2) * raw_z +
    raw_to_voxel(0, 3);
  vy = raw_to_voxel(1, 0) * raw_x +
    raw_to_voxel(1, 1) * raw_y +
    raw_to_voxel(1, 2) * raw_z +
    raw_to_voxel(1, 3);
  vz = raw_to_voxel(2, 0) * raw_x +
    raw_to_voxel(2, 1) * raw_y +
    raw_to_voxel(2, 2) * raw_z +
    raw_to_voxel(2, 3);
}

std::vector<int> voxelized_track_keys(const NumericMatrix &coords,
                                      const NumericMatrix &raw_to_voxel,
                                      int method,
                                      const LogicalVector &axis_flip,
                                      const IntegerVector &vol_dim,
                                      bool interpolate_segments) {
  std::vector<int> keys;
  if (coords.nrow() == 0) return keys;
  if (!interpolate_segments || coords.nrow() == 1) {
    keys.reserve(coords.nrow());
    for (int i = 0; i < coords.nrow(); ++i) {
      int key = -1;
      discretized_key(
        coords(i, 0), coords(i, 1), coords(i, 2), raw_to_voxel,
        method, axis_flip, vol_dim, key
      );
      keys.push_back(key);
    }
    return keys;
  }

  double prev_x;
  double prev_y;
  double prev_z;
  raw_to_voxel_point(coords(0, 0), coords(0, 1), coords(0, 2),
                     raw_to_voxel, prev_x, prev_y, prev_z);
  for (int i = 1; i < coords.nrow(); ++i) {
    double next_x;
    double next_y;
    double next_z;
    raw_to_voxel_point(coords(i, 0), coords(i, 1), coords(i, 2),
                       raw_to_voxel, next_x, next_y, next_z);
    if (std::isfinite(prev_x) && std::isfinite(prev_y) && std::isfinite(prev_z) &&
        std::isfinite(next_x) && std::isfinite(next_y) && std::isfinite(next_z)) {
      const double dx = next_x - prev_x;
      const double dy = next_y - prev_y;
      const double dz = next_z - prev_z;
      const double max_axis_delta = std::max(std::fabs(dx),
        std::max(std::fabs(dy), std::fabs(dz)));
      const int steps = std::max(1, static_cast<int>(std::ceil(max_axis_delta / 0.5)));
      keys.reserve(keys.size() + steps);
      for (int step = 0; step < steps; ++step) {
        const double fraction = static_cast<double>(step) / static_cast<double>(steps);
        int key = -1;
        voxel_key(
          prev_x + fraction * dx,
          prev_y + fraction * dy,
          prev_z + fraction * dz,
          method, axis_flip, vol_dim, key
        );
        keys.push_back(key);
      }
    } else {
      keys.push_back(-1);
    }
    prev_x = next_x;
    prev_y = next_y;
    prev_z = next_z;
  }
  int final_key = -1;
  voxel_key(prev_x, prev_y, prev_z, method, axis_flip, vol_dim, final_key);
  keys.push_back(final_key);
  return keys;
}

enum TargetIntersectionMode {
  TARGET_FULL_STREAMLINE,
  TARGET_ENDPOINT_VOXELS,
  TARGET_LOCAL_SEARCH
};

TargetIntersectionMode target_intersection_mode_code(
    const std::string &mode,
    const std::string &context) {
  if (mode == "full_streamline") return TARGET_FULL_STREAMLINE;
  if (mode == "endpoint_voxels") return TARGET_ENDPOINT_VOXELS;
  if (mode == "local_search") return TARGET_LOCAL_SEARCH;
  stop("Unsupported target_intersection_mode in " + context + ": " + mode);
}

void validate_local_search_parameters(
    double radius_mm,
    const NumericMatrix &voxel_to_world,
    const NumericVector &voxel_radius,
    const std::string &context) {
  if (!std::isfinite(radius_mm) || radius_mm <= 0.0) {
    stop("local_search_radius_mm must be finite and positive in " +
         context + ".");
  }
  if (voxel_to_world.nrow() != 3 || voxel_to_world.ncol() != 3) {
    stop("local_search_voxel_to_world must be a 3-by-3 matrix in " +
         context + ".");
  }
  if (voxel_radius.size() != 3) {
    stop("local_search_voxel_radius must have length 3 in " + context + ".");
  }
  double normalized_basis[3][3];
  for (int axis = 0; axis < 3; ++axis) {
    if (!std::isfinite(voxel_radius[axis]) || voxel_radius[axis] < 0.0) {
      stop("local_search_voxel_radius values must be finite and "
           "non-negative in " + context + ".");
    }
    double row_norm_squared = 0.0;
    for (int component = 0; component < 3; ++component) {
      const double value = voxel_to_world(axis, component);
      if (!std::isfinite(value)) {
        stop("local_search_voxel_to_world values must be finite in " +
             context + ".");
      }
      row_norm_squared += value * value;
    }
    if (!(row_norm_squared > 0.0) || !std::isfinite(row_norm_squared)) {
      stop("local_search_voxel_to_world rows must be finite, non-zero "
           "voxel basis vectors in " + context + ".");
    }
    const double row_norm = std::sqrt(row_norm_squared);
    for (int component = 0; component < 3; ++component) {
      normalized_basis[axis][component] =
        voxel_to_world(axis, component) / row_norm;
    }
  }
  const double normalized_determinant =
    normalized_basis[0][0] *
      (normalized_basis[1][1] * normalized_basis[2][2] -
       normalized_basis[1][2] * normalized_basis[2][1]) -
    normalized_basis[0][1] *
      (normalized_basis[1][0] * normalized_basis[2][2] -
       normalized_basis[1][2] * normalized_basis[2][0]) +
    normalized_basis[0][2] *
      (normalized_basis[1][0] * normalized_basis[2][1] -
       normalized_basis[1][1] * normalized_basis[2][0]);
  if (!std::isfinite(normalized_determinant) ||
      std::fabs(normalized_determinant) <=
        64.0 * std::numeric_limits<double>::epsilon()) {
    stop("local_search_voxel_to_world must contain three non-singular "
         "voxel basis vectors in " + context + ".");
  }
}

int voxel_node_label(int key,
                     const LogicalVector &lut_f,
                     const LogicalMatrix &target_lut) {
  if (key < 0) return -1;
  int label = -1;
  if (lut_f[key] == TRUE) label = 0;
  for (int target = 0; target < target_lut.ncol(); ++target) {
    if (target_lut(key, target) != TRUE) continue;
    if (label >= 0) return -2;
    label = target + 1;
  }
  return label;
}

void endpoint_pair_target_hits(int start_label,
                               int end_label,
                               bool two_endpoints,
                               std::vector<int> &target_hit) {
  if (!two_endpoints || start_label < 0 || end_label < 0) return;
  if (start_label == 0 && end_label > 0 &&
      end_label <= static_cast<int>(target_hit.size())) {
    target_hit[end_label - 1] = 1;
  } else if (end_label == 0 && start_label > 0 &&
             start_label <= static_cast<int>(target_hit.size())) {
    target_hit[start_label - 1] = 1;
  }
}

bool continuous_voxel_point(const NumericMatrix &coords,
                            int row,
                            const NumericMatrix &raw_to_voxel,
                            const LogicalVector &axis_flip,
                            const IntegerVector &vol_dim,
                            double voxel[3]) {
  raw_to_voxel_point(
    coords(row, 0), coords(row, 1), coords(row, 2),
    raw_to_voxel, voxel[0], voxel[1], voxel[2]
  );
  for (int axis = 0; axis < 3; ++axis) {
    if (!std::isfinite(voxel[axis])) return false;
    if (axis_flip[axis] == TRUE) {
      voxel[axis] =
        static_cast<double>(vol_dim[axis]) - voxel[axis] + 1.0;
      if (!std::isfinite(voxel[axis])) return false;
    }
  }
  return true;
}

double squared_distance_tolerance(double left, double right) {
  const double scale = std::max(
    1.0, std::max(std::fabs(left), std::fabs(right))
  );
  return 1e-12 * scale;
}

double local_node_distance_squared(
    int x,
    int y,
    int z,
    const double endpoint_voxel[3],
    const NumericMatrix &voxel_to_world) {
  const double delta[3] = {
    static_cast<double>(x) - endpoint_voxel[0],
    static_cast<double>(y) - endpoint_voxel[1],
    static_cast<double>(z) - endpoint_voxel[2]
  };
  double scanner_delta[3] = {0.0, 0.0, 0.0};
  for (int axis = 0; axis < 3; ++axis) {
    for (int component = 0; component < 3; ++component) {
      scanner_delta[component] +=
        delta[axis] * voxel_to_world(axis, component);
    }
  }
  return scanner_delta[0] * scanner_delta[0] +
    scanner_delta[1] * scanner_delta[1] +
    scanner_delta[2] * scanner_delta[2];
}

double scanner_metric_lower_bound(
    const NumericMatrix &voxel_to_world) {
  // If singular values are s1 >= s2 >= s3, then
  // s1 * s2 <= ||B||_F^2 / 2 and therefore
  // s3^2 >= (2 * |det(B)| / ||B||_F^2)^2.  Normalize first and
  // widen determinant/Frobenius rounding conservatively so this remains
  // suitable for pruning reflected, sheared, and oblique bases.
  long double scale = 0.0L;
  for (int axis = 0; axis < 3; ++axis) {
    for (int component = 0; component < 3; ++component) {
      scale = std::max(
        scale,
        std::fabs(
          static_cast<long double>(voxel_to_world(axis, component))
        )
      );
    }
  }
  if (!(scale > 0.0L) || !std::isfinite(scale)) return 0.0;

  long double basis[3][3];
  long double frobenius_squared = 0.0L;
  for (int axis = 0; axis < 3; ++axis) {
    for (int component = 0; component < 3; ++component) {
      basis[axis][component] =
        static_cast<long double>(voxel_to_world(axis, component)) / scale;
      frobenius_squared +=
        basis[axis][component] * basis[axis][component];
    }
  }
  if (!(frobenius_squared > 0.0L) ||
      !std::isfinite(frobenius_squared)) {
    return 0.0;
  }

  const long double positive_terms[3] = {
    basis[0][0] * basis[1][1] * basis[2][2],
    basis[0][1] * basis[1][2] * basis[2][0],
    basis[0][2] * basis[1][0] * basis[2][1]
  };
  const long double negative_terms[3] = {
    basis[0][2] * basis[1][1] * basis[2][0],
    basis[0][1] * basis[1][0] * basis[2][2],
    basis[0][0] * basis[1][2] * basis[2][1]
  };
  long double determinant = 0.0L;
  long double determinant_term_sum = 0.0L;
  for (int term = 0; term < 3; ++term) {
    determinant += positive_terms[term] - negative_terms[term];
    determinant_term_sum +=
      std::fabs(positive_terms[term]) +
      std::fabs(negative_terms[term]);
  }

  const long double calculation_error =
    64.0L * std::numeric_limits<long double>::epsilon();
  const long double safe_abs_determinant = std::max(
    0.0L,
    std::fabs(determinant) -
      calculation_error * std::max(1.0L, determinant_term_sum)
  );
  if (!(safe_abs_determinant > 0.0L)) return 0.0;

  const long double frobenius_upper =
    frobenius_squared * (1.0L + calculation_error);
  const long double normalized_lower =
    4.0L * safe_abs_determinant * safe_abs_determinant /
    (frobenius_upper * frobenius_upper);
  const long double scaled_lower =
    normalized_lower * scale * scale;
  if (!(scaled_lower > 0.0L) ||
      !std::isfinite(scaled_lower) ||
      scaled_lower >
        static_cast<long double>(std::numeric_limits<double>::max())) {
    return 0.0;
  }

  return std::nextafter(
    static_cast<double>(scaled_lower), 0.0
  );
}

double squared_tolerance_relevance_limit(double reference) {
  // Any d satisfying d <= reference +
  // 1e-12 * max(1, d, reference) lies below this widened limit.
  if (!std::isfinite(reference)) return R_PosInf;
  const double scale = std::max(1.0, std::fabs(reference));
  const double limit = reference + 2e-12 * scale;
  if (!std::isfinite(limit)) return R_PosInf;
  return std::nextafter(limit, R_PosInf);
}

int nearest_local_node_label(
    const double endpoint_voxel[3],
    double radius_mm,
    const NumericMatrix &voxel_to_world,
    const NumericVector &voxel_radius,
    const IntegerVector &vol_dim,
    const LogicalVector &lut_f,
    const LogicalMatrix &target_lut,
    double *distance_mm = NULL,
    bool *ambiguous = NULL) {
  if (distance_mm != NULL) *distance_mm = NA_REAL;
  if (ambiguous != NULL) *ambiguous = false;
  int lower[3];
  int upper[3];
  for (int axis = 0; axis < 3; ++axis) {
    const double candidate_lower =
      std::ceil(endpoint_voxel[axis] - voxel_radius[axis]);
    const double candidate_upper =
      std::floor(endpoint_voxel[axis] + voxel_radius[axis]);
    const double clipped_lower = std::max(1.0, candidate_lower);
    const double clipped_upper = std::min(
      static_cast<double>(vol_dim[axis]), candidate_upper
    );
    if (!std::isfinite(clipped_lower) ||
        !std::isfinite(clipped_upper) ||
        clipped_lower > clipped_upper) {
      return -1;
    }
    lower[axis] = static_cast<int>(clipped_lower);
    upper[axis] = static_cast<int>(clipped_upper);
  }

  const double radius_squared = radius_mm * radius_mm;
  const double metric_lower_bound =
    scanner_metric_lower_bound(voxel_to_world);

  // Phase one visits each lattice point once in centre-out Chebyshev shells.
  // It discovers a tight in-radius upper bound without changing the final
  // reduction order used for near-distance ties.
  int center[3];
  double center_offset_squared = 0.0;
  int max_shell = 0;
  for (int axis = 0; axis < 3; ++axis) {
    const double rounded = std::nearbyint(endpoint_voxel[axis]);
    if (rounded <= static_cast<double>(lower[axis])) {
      center[axis] = lower[axis];
    } else if (rounded >= static_cast<double>(upper[axis])) {
      center[axis] = upper[axis];
    } else {
      center[axis] = static_cast<int>(rounded);
    }
    const double offset =
      endpoint_voxel[axis] - static_cast<double>(center[axis]);
    center_offset_squared += offset * offset;
    max_shell = std::max(
      max_shell,
      std::max(
        center[axis] - lower[axis],
        upper[axis] - center[axis]
      )
    );
  }
  const double center_offset = std::sqrt(center_offset_squared);

  double provisional_distance_squared = R_PosInf;
  auto inspect_provisional = [&](int x, int y, int z) {
    const int key = encode_key(x, y, z, vol_dim);
    if (voxel_node_label(key, lut_f, target_lut) == -1) return;
    const double distance_squared = local_node_distance_squared(
      x, y, z, endpoint_voxel, voxel_to_world
    );
    if (!std::isfinite(distance_squared)) return;
    const double radius_tolerance =
      squared_distance_tolerance(distance_squared, radius_squared);
    if (distance_squared > radius_squared + radius_tolerance) return;
    if (distance_squared < provisional_distance_squared) {
      provisional_distance_squared = distance_squared;
    }
  };

  for (int shell = 0; shell <= max_shell; ++shell) {
    const int x_low = std::max(lower[0], center[0] - shell);
    const int x_high = std::min(upper[0], center[0] + shell);
    if (shell == 0) {
      inspect_provisional(center[0], center[1], center[2]);
    } else {
      const int z_faces[2] = {
        center[2] - shell, center[2] + shell
      };
      for (int face = 0; face < 2; ++face) {
        const int z = z_faces[face];
        if (z < lower[2] || z > upper[2]) continue;
        const int y_low = std::max(lower[1], center[1] - shell);
        const int y_high = std::min(upper[1], center[1] + shell);
        for (int y = y_low; y <= y_high; ++y) {
          for (int x = x_low; x <= x_high; ++x) {
            inspect_provisional(x, y, z);
          }
        }
      }

      const int z_inner_low =
        std::max(lower[2], center[2] - shell + 1);
      const int z_inner_high =
        std::min(upper[2], center[2] + shell - 1);
      const int y_faces[2] = {
        center[1] - shell, center[1] + shell
      };
      for (int z = z_inner_low; z <= z_inner_high; ++z) {
        for (int face = 0; face < 2; ++face) {
          const int y = y_faces[face];
          if (y < lower[1] || y > upper[1]) continue;
          for (int x = x_low; x <= x_high; ++x) {
            inspect_provisional(x, y, z);
          }
        }
      }

      const int y_inner_low =
        std::max(lower[1], center[1] - shell + 1);
      const int y_inner_high =
        std::min(upper[1], center[1] + shell - 1);
      const int x_faces[2] = {
        center[0] - shell, center[0] + shell
      };
      for (int z = z_inner_low; z <= z_inner_high; ++z) {
        for (int y = y_inner_low; y <= y_inner_high; ++y) {
          for (int face = 0; face < 2; ++face) {
            const int x = x_faces[face];
            if (x < lower[0] || x > upper[0]) continue;
            inspect_provisional(x, y, z);
          }
        }
      }
    }

    if (shell == max_shell) break;
    const double lower_bound_scale = std::max(
      1.0,
      std::max(
        static_cast<double>(shell + 1),
        center_offset
      )
    );
    const double index_distance_lower = std::max(
      0.0,
      static_cast<double>(shell + 1) - center_offset -
        64.0 * std::numeric_limits<double>::epsilon() *
          lower_bound_scale
    );
    double remaining_distance_lower =
      metric_lower_bound *
      index_distance_lower * index_distance_lower;
    if (remaining_distance_lower > 0.0 &&
        std::isfinite(remaining_distance_lower)) {
      remaining_distance_lower =
        std::nextafter(remaining_distance_lower, 0.0);
    } else {
      remaining_distance_lower = 0.0;
    }

    if (std::isfinite(provisional_distance_squared)) {
      const double nearest_tolerance = squared_distance_tolerance(
        remaining_distance_lower, provisional_distance_squared
      );
      if (remaining_distance_lower >
          provisional_distance_squared + nearest_tolerance) {
        break;
      }
    } else {
      const double radius_tolerance = squared_distance_tolerance(
        remaining_distance_lower, radius_squared
      );
      if (remaining_distance_lower > radius_squared + radius_tolerance) {
        return -1;
      }
    }
  }
  if (!std::isfinite(provisional_distance_squared)) return -1;

  const double relevance_limit = squared_tolerance_relevance_limit(
    provisional_distance_squared
  );

  // Phase two performs the original z/y/x lexicographic reduction.  Since
  // d' G d >= metric_lower_bound * ||d||^2, every relevant lattice point is
  // inside the index-space ball below.  One extra lattice voxel on each face
  // conservatively absorbs floating-point projection error.  If the metric
  // is too ill-conditioned for a positive certified bound, retain the full
  // original box and sacrifice only pruning, not semantics.
  int exact_lower[3];
  int exact_upper[3];
  const double exact_extent =
    metric_lower_bound > 0.0 && std::isfinite(relevance_limit) ?
      std::sqrt(relevance_limit / metric_lower_bound) :
      R_PosInf;
  for (int axis = 0; axis < 3; ++axis) {
    if (!std::isfinite(exact_extent)) {
      exact_lower[axis] = lower[axis];
      exact_upper[axis] = upper[axis];
      continue;
    }
    const double candidate_lower =
      std::ceil(endpoint_voxel[axis] - exact_extent) - 1.0;
    const double candidate_upper =
      std::floor(endpoint_voxel[axis] + exact_extent) + 1.0;
    exact_lower[axis] = std::max(
      lower[axis],
      static_cast<int>(std::max(
        static_cast<double>(lower[axis]), candidate_lower
      ))
    );
    exact_upper[axis] = std::min(
      upper[axis],
      static_cast<int>(std::min(
        static_cast<double>(upper[axis]), candidate_upper
      ))
    );
  }

  double best_distance_squared = R_PosInf;
  int best_label = -1;
  bool best_ambiguous = false;
  for (int z = exact_lower[2]; z <= exact_upper[2]; ++z) {
    for (int y = exact_lower[1]; y <= exact_upper[1]; ++y) {
      for (int x = exact_lower[0]; x <= exact_upper[0]; ++x) {
        const int key = encode_key(x, y, z, vol_dim);
        const int label = voxel_node_label(key, lut_f, target_lut);
        if (label == -1) continue;

        const double distance_squared = local_node_distance_squared(
          x, y, z, endpoint_voxel, voxel_to_world
        );
        if (!std::isfinite(distance_squared)) continue;
        const double radius_tolerance =
          squared_distance_tolerance(distance_squared, radius_squared);
        if (distance_squared > radius_squared + radius_tolerance) continue;

        if (!std::isfinite(best_distance_squared)) {
          best_distance_squared = distance_squared;
          best_label = label;
          best_ambiguous = label == -2;
          continue;
        }
        const double nearest_tolerance = squared_distance_tolerance(
          distance_squared, best_distance_squared
        );
        if (distance_squared < best_distance_squared - nearest_tolerance) {
          best_distance_squared = distance_squared;
          best_label = label;
          best_ambiguous = label == -2;
        } else if (std::fabs(distance_squared - best_distance_squared) <=
                   nearest_tolerance) {
          if (label == -2 || best_label == -2 || label != best_label) {
            best_ambiguous = true;
          }
        }
      }
    }
  }
  if (distance_mm != NULL && std::isfinite(best_distance_squared)) {
    *distance_mm = std::sqrt(best_distance_squared);
  }
  if (best_ambiguous) {
    if (ambiguous != NULL) *ambiguous = true;
    return -2;
  }
  if (best_label < 0) return -1;
  return best_label;
}

std::vector<int> primary_target_hits(
    TargetIntersectionMode mode,
    const std::vector<int> &full_path_target_hit,
    const std::vector<int> &chosen_keys,
    const NumericMatrix &coords,
    const NumericMatrix &raw_to_voxel,
    const LogicalVector &axis_flip,
    const IntegerVector &vol_dim,
    const LogicalVector &lut_f,
    const LogicalMatrix &target_lut,
    double local_search_radius_mm,
    const NumericMatrix &local_search_voxel_to_world,
    const NumericVector &local_search_voxel_radius,
    int *start_label_output = NULL,
    int *end_label_output = NULL,
    double *start_distance_output = NULL,
    double *end_distance_output = NULL,
    bool *start_ambiguous_output = NULL,
    bool *end_ambiguous_output = NULL) {
  if (start_label_output != NULL) *start_label_output = -3;
  if (end_label_output != NULL) *end_label_output = -3;
  if (start_distance_output != NULL) *start_distance_output = NA_REAL;
  if (end_distance_output != NULL) *end_distance_output = NA_REAL;
  if (start_ambiguous_output != NULL) *start_ambiguous_output = false;
  if (end_ambiguous_output != NULL) *end_ambiguous_output = false;
  if (mode == TARGET_FULL_STREAMLINE) return full_path_target_hit;

  std::vector<int> target_hit(target_lut.ncol(), 0);
  const bool two_endpoints = coords.nrow() > 1;
  if (!two_endpoints) return target_hit;

  int start_label = -1;
  int end_label = -1;
  if (mode == TARGET_ENDPOINT_VOXELS) {
    start_label = voxel_node_label(chosen_keys.front(), lut_f, target_lut);
    end_label = voxel_node_label(chosen_keys.back(), lut_f, target_lut);
  } else {
    double start_voxel[3];
    double end_voxel[3];
    if (continuous_voxel_point(
          coords, 0, raw_to_voxel, axis_flip, vol_dim, start_voxel)) {
      start_label = nearest_local_node_label(
        start_voxel, local_search_radius_mm,
        local_search_voxel_to_world, local_search_voxel_radius,
        vol_dim, lut_f, target_lut,
        start_distance_output, start_ambiguous_output
      );
    }
    if (continuous_voxel_point(
          coords, coords.nrow() - 1, raw_to_voxel,
          axis_flip, vol_dim, end_voxel)) {
      end_label = nearest_local_node_label(
        end_voxel, local_search_radius_mm,
        local_search_voxel_to_world, local_search_voxel_radius,
        vol_dim, lut_f, target_lut,
        end_distance_output, end_ambiguous_output
      );
    }
  }
  if (start_label_output != NULL) *start_label_output = start_label;
  if (end_label_output != NULL) *end_label_output = end_label;
  endpoint_pair_target_hits(
    start_label, end_label, two_endpoints, target_hit
  );
  return target_hit;
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
List ftq_density_chunk_cpp(List tracks,
                           NumericMatrix raw_to_voxel,
                           std::string method,
                           LogicalVector axis_flip,
                           IntegerVector vol_dim,
                           std::string streamline_voxelization_mode) {
  if (raw_to_voxel.nrow() < 3 || raw_to_voxel.ncol() < 4) {
    stop("raw_to_voxel must have at least 3 rows and 4 columns.");
  }
  if (axis_flip.size() != 3 || vol_dim.size() < 3) {
    stop("axis_flip and vol_dim must have length 3.");
  }
  const std::vector<int> codes = method_codes(CharacterVector::create(method));
  const int method_code = codes[0];
  const bool interpolate_segments = streamline_voxelization_mode == "segment_interpolation";
  if (streamline_voxelization_mode != "sample_points" && !interpolate_segments) {
    stop("Unsupported streamline_voxelization_mode in C++ density builder: " +
         streamline_voxelization_mode);
  }

  const int nx = vol_dim[0];
  const int ny = vol_dim[1];
  const int nz = vol_dim[2];
  IntegerVector density(static_cast<R_xlen_t>(nx) * ny * nz);
  double stored_points = 0.0;
  double voxelization_points = 0.0;
  double out_of_bounds_points = 0.0;
  double occupancy_events = 0.0;
  int min_x = nx + 1;
  int min_y = ny + 1;
  int min_z = nz + 1;
  int max_x = 0;
  int max_y = 0;
  int max_z = 0;

  for (R_xlen_t t = 0; t < tracks.size(); ++t) {
    NumericMatrix coords = track_coords(tracks[t]);
    if (coords.ncol() < 3 || coords.nrow() == 0) continue;
    stored_points += coords.nrow();
    std::vector<int> keys = voxelized_track_keys(
      coords, raw_to_voxel, method_code, axis_flip, vol_dim,
      interpolate_segments
    );
    voxelization_points += keys.size();
    std::unordered_set<int> seen;
    seen.reserve(keys.size());
    for (std::size_t i = 0; i < keys.size(); ++i) {
      const int key = keys[i];
      if (key < 0) {
        out_of_bounds_points += 1.0;
        continue;
      }
      if (!seen.insert(key).second) continue;
      density[key] += 1;
      occupancy_events += 1.0;
      const int x = key % nx + 1;
      const int y = (key / nx) % ny + 1;
      const int z = key / (nx * ny) + 1;
      if (x < min_x) min_x = x;
      if (y < min_y) min_y = y;
      if (z < min_z) min_z = z;
      if (x > max_x) max_x = x;
      if (y > max_y) max_y = y;
      if (z > max_z) max_z = z;
    }
  }

  return List::create(
    _["density"] = density,
    _["processed"] = static_cast<int>(tracks.size()),
    _["stored_points"] = stored_points,
    _["voxelization_points"] = voxelization_points,
    _["out_of_bounds_points"] = out_of_bounds_points,
    _["occupancy_events"] = occupancy_events,
    _["voxel_min"] = IntegerVector::create(min_x, min_y, min_z),
    _["voxel_max"] = IntegerVector::create(max_x, max_y, max_z)
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

// [[Rcpp::export]]
List ftq_count_chunk_cpp(List tracks,
                         NumericMatrix raw_to_voxel,
                         CharacterVector methods,
                         LogicalVector axis_flip,
                         IntegerVector vol_dim,
                         IntegerVector f_index,
                         LogicalVector lut_f,
                         LogicalVector lut_local,
                         bool has_local,
                         LogicalMatrix target_lut,
                         int n_fnode,
                         std::string focal_assignment_mode,
                         std::string streamline_voxelization_mode,
                         std::string target_intersection_mode,
                         double local_search_radius_mm,
                         NumericMatrix local_search_voxel_to_world,
                         NumericVector local_search_voxel_radius,
                         bool endpoint_analysis,
                         LogicalVector lut_f_dilated,
                         LogicalMatrix target_lut_dilated,
                         LogicalMatrix context_lut) {
  if (raw_to_voxel.nrow() < 3 || raw_to_voxel.ncol() < 4) {
    stop("raw_to_voxel must have at least 3 rows and 4 columns.");
  }
  if (axis_flip.size() != 3 || vol_dim.size() < 3) {
    stop("axis_flip and vol_dim must have length 3.");
  }
  const bool all_focal = focal_assignment_mode == "all_focal_voxels";
  const bool first_intersection = focal_assignment_mode == "first_intersection";
	  if (!all_focal && !first_intersection) {
	    stop("Unsupported focal_assignment_mode in C++ backend: " + focal_assignment_mode);
	  }
	  const bool interpolate_segments = streamline_voxelization_mode == "segment_interpolation";
  if (streamline_voxelization_mode != "sample_points" && !interpolate_segments) {
    stop("Unsupported streamline_voxelization_mode in C++ backend: " + streamline_voxelization_mode);
  }
  const TargetIntersectionMode target_mode = target_intersection_mode_code(
    target_intersection_mode, "C++ backend"
  );
  if (target_mode == TARGET_LOCAL_SEARCH) {
    validate_local_search_parameters(
      local_search_radius_mm, local_search_voxel_to_world,
      local_search_voxel_radius, "C++ backend"
    );
  }

  const std::vector<int> m_codes = method_codes(methods);
  if (m_codes.empty()) stop("At least one discretization method is required.");
  const int n_targets = target_lut.ncol();
  const int n_voxels = vol_dim[0] * vol_dim[1] * vol_dim[2];
  const int n_contexts = context_lut.ncol();
  if (target_lut.nrow() != n_voxels || lut_f.size() != n_voxels) {
    stop("Focal and target lookup tables must match the reference grid.");
  }
  if (endpoint_analysis &&
      (lut_f_dilated.size() != n_voxels ||
       target_lut_dilated.nrow() != n_voxels ||
       target_lut_dilated.ncol() != n_targets ||
       context_lut.nrow() != n_voxels)) {
    stop("Endpoint lookup tables must match the reference grid and target count.");
  }

  IntegerVector tracks_counts(n_fnode);
  IntegerVector local_counts(n_fnode);
  IntegerVector leave_no_target_counts(n_fnode);
  IntegerMatrix target_counts(n_fnode, n_targets);
  IntegerVector target_unique_counts(n_targets);
  IntegerVector full_path_target_unique_counts(n_targets);
  std::vector<int> local_search_chunk_streamline_index;
  std::vector<int> local_search_start_label;
  std::vector<int> local_search_end_label;
  std::vector<double> local_search_start_distance;
  std::vector<double> local_search_end_distance;
  std::vector<int> local_search_start_ambiguous;
  std::vector<int> local_search_end_ambiguous;
  std::vector<int> local_search_paired_target_index;
  std::vector<double> local_search_paired_focal_distance;
  std::vector<double> local_search_paired_target_distance;
  std::vector<double> local_search_paired_max_distance;
  std::vector<int> local_search_paired_full_path;
  int processed = 0;
  double unique_focal_streamlines = 0.0;
	  double unique_local_streamlines = 0.0;
	  double unique_leave_no_target_streamlines = 0.0;
	  double unique_with_target_streamlines = 0.0;
	  double unique_local_target_positive_streamlines = 0.0;
	  double unique_local_target_negative_streamlines = 0.0;
	  double unique_nonlocal_target_positive_streamlines = 0.0;
	  double unique_nonlocal_target_negative_streamlines = 0.0;
	  double unique_any_target_positive_streamlines = 0.0;
	  double unique_any_target_negative_streamlines = 0.0;
	  double total_focal_associations = 0.0;
  double total_focal_target_associations = 0.0;

  double endpoint_valid_streamlines = 0.0;
  double endpoint_count = 0.0;
  double endpoint_out_of_bounds = 0.0;
  double endpoint_nonfinite = 0.0;
  double endpoint_outside_one_voxel = 0.0;
  double endpoint_outside_more_than_one_voxel = 0.0;
  NumericVector endpoint_out_of_bounds_low(3);
  NumericVector endpoint_out_of_bounds_high(3);
  NumericMatrix endpoint_out_of_bounds_low_distance(3, 2);
  NumericMatrix endpoint_out_of_bounds_high_distance(3, 2);
  double grid_any_endpoint_inside = 0.0;
  double grid_both_endpoints_inside = 0.0;
  NumericVector context_endpoint_inside(n_contexts);
  NumericVector context_any_endpoint_inside(n_contexts);
  NumericVector context_both_endpoints_inside(n_contexts);
  double focal_endpoint_streamlines_exact = 0.0;
  double focal_endpoint_streamlines_dilated1 = 0.0;
  NumericVector paired_endpoint_streamlines_exact(n_targets);
  NumericVector paired_endpoint_streamlines_dilated1(n_targets);
  NumericVector paired_intersection_streamlines_exact(n_targets);
  NumericVector paired_intersection_streamlines_dilated1(n_targets);
  NumericVector focal_endpoint_target_intersection_exact(n_targets);
  NumericVector focal_endpoint_target_intersection_dilated1(n_targets);
  NumericMatrix endpoint_class_exact(n_targets, 5);
  NumericMatrix endpoint_class_dilated1(n_targets, 5);
  std::vector<std::unordered_map<int, double> > penetration_hist(n_targets);

  const auto record_context_endpoint = [&](const NumericMatrix &coords,
                                           int row,
                                           int &key) {
    double vx;
    double vy;
    double vz;
    raw_to_voxel_point(
      coords(row, 0), coords(row, 1), coords(row, 2),
      raw_to_voxel, vx, vy, vz
    );
    const double values[3] = {vx, vy, vz};
    int discrete[3] = {0, 0, 0};
    bool finite = true;
    int max_distance = 0;
    for (int axis = 0; axis < 3; ++axis) {
      if (!std::isfinite(values[axis])) {
        finite = false;
        continue;
      }
      double value;
      if (m_codes[0] == 0) {
        value = std::nearbyint(values[axis]);
      } else if (m_codes[0] == 1) {
        value = std::floor(values[axis]);
      } else {
        value = std::ceil(values[axis]);
      }
      discrete[axis] = static_cast<int>(value);
      if (axis_flip[axis] == TRUE) {
        discrete[axis] = vol_dim[axis] - discrete[axis] + 1;
      }
      if (discrete[axis] < 1) {
        const int distance = 1 - discrete[axis];
        endpoint_out_of_bounds_low[axis] += 1.0;
        endpoint_out_of_bounds_low_distance(axis, distance == 1 ? 0 : 1) += 1.0;
        max_distance = std::max(max_distance, distance);
      } else if (discrete[axis] > vol_dim[axis]) {
        const int distance = discrete[axis] - vol_dim[axis];
        endpoint_out_of_bounds_high[axis] += 1.0;
        endpoint_out_of_bounds_high_distance(axis, distance == 1 ? 0 : 1) += 1.0;
        max_distance = std::max(max_distance, distance);
      }
    }
    endpoint_count += 1.0;
    if (!finite) {
      key = -1;
      endpoint_nonfinite += 1.0;
      endpoint_out_of_bounds += 1.0;
      return false;
    }
    key = encode_key(discrete[0], discrete[1], discrete[2], vol_dim);
    if (key < 0) {
      endpoint_out_of_bounds += 1.0;
      if (max_distance == 1) {
        endpoint_outside_one_voxel += 1.0;
      } else if (max_distance > 1) {
        endpoint_outside_more_than_one_voxel += 1.0;
      }
      return false;
    }
    return true;
  };

  const auto classify_endpoint_pair = [](bool start_focal,
                                         bool end_focal,
                                         bool start_target,
                                         bool end_target,
                                         bool two_endpoints,
                                         bool &paired,
                                         bool &focal_any,
                                         bool &target_any) {
    focal_any = start_focal || end_focal;
    target_any = start_target || end_target;
    const bool shared =
      (start_focal && start_target) || (end_focal && end_target);
    paired = two_endpoints &&
      ((start_focal && end_target) || (end_focal && start_target));
    if (shared) return 0;
    if (paired) return 1;
    if (target_any) return 2;
    if (focal_any) return 3;
    return 4;
  };

  for (R_xlen_t t = 0; t < tracks.size(); ++t) {
    NumericMatrix coords = track_coords(tracks[t]);
    ++processed;
    if (coords.ncol() < 3 || coords.nrow() == 0) continue;

    const bool two_endpoints = coords.nrow() > 1;
    if (endpoint_analysis) {
      endpoint_valid_streamlines += 1.0;
      int context_start_key = -1;
      int context_end_key = -1;
      const bool start_inside_grid = record_context_endpoint(coords, 0, context_start_key);
      const bool end_inside_grid = two_endpoints ?
        record_context_endpoint(coords, coords.nrow() - 1, context_end_key) :
        start_inside_grid;
      if (!two_endpoints) context_end_key = context_start_key;
      grid_any_endpoint_inside +=
        (start_inside_grid || end_inside_grid) ? 1.0 : 0.0;
      grid_both_endpoints_inside +=
        (start_inside_grid && end_inside_grid) ? 1.0 : 0.0;
      for (int context = 0; context < n_contexts; ++context) {
        const bool start_inside = context_start_key >= 0 &&
          context_lut(context_start_key, context) == TRUE;
        const bool end_inside = context_end_key >= 0 &&
          context_lut(context_end_key, context) == TRUE;
        context_endpoint_inside[context] += start_inside ? 1.0 : 0.0;
        if (two_endpoints) context_endpoint_inside[context] += end_inside ? 1.0 : 0.0;
        context_any_endpoint_inside[context] +=
          (start_inside || end_inside) ? 1.0 : 0.0;
        context_both_endpoints_inside[context] +=
          (start_inside && end_inside) ? 1.0 : 0.0;
      }
    }

    std::vector<int> chosen_keys;
    std::vector<int> focal_indices;

    for (std::size_t mi = 0; mi < m_codes.size(); ++mi) {
      std::vector<int> keys = voxelized_track_keys(
        coords, raw_to_voxel, m_codes[mi], axis_flip, vol_dim,
        interpolate_segments
      );
      std::vector<int> method_focal_indices;
      std::unordered_set<int> seen_focal_indices;
      int first_hit_origin = -1;

      for (std::size_t i = 0; i < keys.size(); ++i) {
        const int key = keys[i];
        if (key >= 0 && lut_f[key] == TRUE) {
          const int idx = f_index[key];
          if (idx == NA_INTEGER) continue;
          const int row_idx = idx - 1;
          if (row_idx < 0 || row_idx >= n_fnode) continue;
          if (first_intersection) {
            if (first_hit_origin < 0) first_hit_origin = row_idx;
          } else if (seen_focal_indices.insert(row_idx).second) {
            method_focal_indices.push_back(row_idx);
          }
        }
      }

      if (first_intersection && first_hit_origin >= 0) {
        method_focal_indices.push_back(first_hit_origin);
      }
      if (!method_focal_indices.empty()) {
        focal_indices.swap(method_focal_indices);
        chosen_keys.swap(keys);
        break;
      }
    }

    if (focal_indices.empty()) continue;

    int local_flag = 0;
    if (has_local) {
      bool all_valid = true;
      bool all_local = true;
      for (std::size_t i = 0; i < chosen_keys.size(); ++i) {
        const int key = chosen_keys[i];
        if (key < 0) {
          all_valid = false;
          all_local = false;
          break;
        }
        if (lut_local[key] != TRUE) all_local = false;
      }
      if (all_valid && all_local) local_flag = 1;
    }

    std::vector<int> full_path_target_hit(n_targets, 0);
    std::vector<int> target_voxel_counts(n_targets, 0);
    if (n_targets > 0) {
      std::unordered_set<int> seen;
      seen.reserve(chosen_keys.size());
      for (std::size_t i = 0; i < chosen_keys.size(); ++i) {
        const int key = chosen_keys[i];
        if (key < 0) continue;
        if (!seen.insert(key).second) continue;
        for (int j = 0; j < n_targets; ++j) {
          if (target_lut(key, j) == TRUE) {
            full_path_target_hit[j] = 1;
            if (endpoint_analysis) target_voxel_counts[j] += 1;
          }
        }
      }
    }
    int local_search_start_code = -3;
    int local_search_end_code = -3;
    double local_search_start_mm = NA_REAL;
    double local_search_end_mm = NA_REAL;
    bool local_search_start_is_ambiguous = false;
    bool local_search_end_is_ambiguous = false;
    std::vector<int> target_hit = primary_target_hits(
      target_mode, full_path_target_hit, chosen_keys, coords,
      raw_to_voxel, axis_flip, vol_dim, lut_f, target_lut,
      local_search_radius_mm, local_search_voxel_to_world,
      local_search_voxel_radius,
      &local_search_start_code, &local_search_end_code,
      &local_search_start_mm, &local_search_end_mm,
      &local_search_start_is_ambiguous,
      &local_search_end_is_ambiguous
    );

    if (target_mode == TARGET_LOCAL_SEARCH) {
      int paired_target_index = NA_INTEGER;
      int paired_full_path = NA_LOGICAL;
      double paired_focal_distance = NA_REAL;
      double paired_target_distance = NA_REAL;
      double paired_max_distance = NA_REAL;
      for (int target = 0; target < n_targets; ++target) {
        if (!target_hit[target]) continue;
        paired_target_index = target + 1;
        paired_full_path = full_path_target_hit[target] ? TRUE : FALSE;
        if (local_search_start_code == 0 &&
            local_search_end_code == target + 1) {
          paired_focal_distance = local_search_start_mm;
          paired_target_distance = local_search_end_mm;
        } else if (local_search_end_code == 0 &&
                   local_search_start_code == target + 1) {
          paired_focal_distance = local_search_end_mm;
          paired_target_distance = local_search_start_mm;
        }
        if (std::isfinite(paired_focal_distance) &&
            std::isfinite(paired_target_distance)) {
          paired_max_distance = std::max(
            paired_focal_distance,
            paired_target_distance
          );
        }
        break;
      }
      local_search_chunk_streamline_index.push_back(
        static_cast<int>(t) + 1
      );
      local_search_start_label.push_back(local_search_start_code);
      local_search_end_label.push_back(local_search_end_code);
      local_search_start_distance.push_back(local_search_start_mm);
      local_search_end_distance.push_back(local_search_end_mm);
      local_search_start_ambiguous.push_back(
        local_search_start_is_ambiguous ? TRUE : FALSE
      );
      local_search_end_ambiguous.push_back(
        local_search_end_is_ambiguous ? TRUE : FALSE
      );
      local_search_paired_target_index.push_back(paired_target_index);
      local_search_paired_focal_distance.push_back(paired_focal_distance);
      local_search_paired_target_distance.push_back(paired_target_distance);
      local_search_paired_max_distance.push_back(paired_max_distance);
      local_search_paired_full_path.push_back(paired_full_path);
    }

    if (endpoint_analysis) {
      const int start_key = chosen_keys.front();
      const int end_key = two_endpoints ? chosen_keys.back() : start_key;
      const bool start_focal_exact = start_key >= 0 && lut_f[start_key] == TRUE;
      const bool end_focal_exact = end_key >= 0 && lut_f[end_key] == TRUE;
      const bool start_focal_d1 = start_key >= 0 && lut_f_dilated[start_key] == TRUE;
      const bool end_focal_d1 = end_key >= 0 && lut_f_dilated[end_key] == TRUE;
      const bool focal_any_exact = start_focal_exact || end_focal_exact;
      const bool focal_any_d1 = start_focal_d1 || end_focal_d1;
      focal_endpoint_streamlines_exact += focal_any_exact ? 1.0 : 0.0;
      focal_endpoint_streamlines_dilated1 += focal_any_d1 ? 1.0 : 0.0;

      for (int j = 0; j < n_targets; ++j) {
        const bool start_target_exact =
          start_key >= 0 && target_lut(start_key, j) == TRUE;
        const bool end_target_exact =
          end_key >= 0 && target_lut(end_key, j) == TRUE;
        const bool start_target_d1 =
          start_key >= 0 && target_lut_dilated(start_key, j) == TRUE;
        const bool end_target_d1 =
          end_key >= 0 && target_lut_dilated(end_key, j) == TRUE;
        bool paired_exact = false;
        bool endpoint_focal_exact = false;
        bool endpoint_target_exact = false;
        bool paired_d1 = false;
        bool endpoint_focal_d1 = false;
        bool endpoint_target_d1 = false;
        const int class_exact = classify_endpoint_pair(
          start_focal_exact, end_focal_exact,
          start_target_exact, end_target_exact, two_endpoints,
          paired_exact, endpoint_focal_exact, endpoint_target_exact
        );
        const int class_d1 = classify_endpoint_pair(
          start_focal_d1, end_focal_d1,
          start_target_d1, end_target_d1, two_endpoints,
          paired_d1, endpoint_focal_d1, endpoint_target_d1
        );
        paired_endpoint_streamlines_exact[j] += paired_exact ? 1.0 : 0.0;
        paired_endpoint_streamlines_dilated1[j] += paired_d1 ? 1.0 : 0.0;
        if (!full_path_target_hit[j]) continue;
        endpoint_class_exact(j, class_exact) += 1.0;
        endpoint_class_dilated1(j, class_d1) += 1.0;
        paired_intersection_streamlines_exact[j] += paired_exact ? 1.0 : 0.0;
        paired_intersection_streamlines_dilated1[j] += paired_d1 ? 1.0 : 0.0;
        focal_endpoint_target_intersection_exact[j] +=
          endpoint_focal_exact ? 1.0 : 0.0;
        focal_endpoint_target_intersection_dilated1[j] +=
          endpoint_focal_d1 ? 1.0 : 0.0;
        penetration_hist[j][target_voxel_counts[j]] += 1.0;
      }
    }

    bool any_target = false;
    for (int j = 0; j < n_targets; ++j) {
      if (full_path_target_hit[j]) {
        full_path_target_unique_counts[j] += 1;
      }
      if (target_hit[j]) {
        any_target = true;
        target_unique_counts[j] += 1;
      }
    }
    unique_focal_streamlines += 1.0;
    unique_local_streamlines += static_cast<double>(local_flag);
	    if (any_target) {
	      unique_with_target_streamlines += 1.0;
	    } else if (local_flag == 0) {
	      unique_leave_no_target_streamlines += 1.0;
	    }
	    if (local_flag == 1 && any_target) {
	      unique_local_target_positive_streamlines += 1.0;
	    } else if (local_flag == 1 && !any_target) {
	      unique_local_target_negative_streamlines += 1.0;
	    } else if (local_flag == 0 && any_target) {
	      unique_nonlocal_target_positive_streamlines += 1.0;
	    } else {
	      unique_nonlocal_target_negative_streamlines += 1.0;
	    }
	    unique_any_target_positive_streamlines += any_target ? 1.0 : 0.0;
	    unique_any_target_negative_streamlines += any_target ? 0.0 : 1.0;
	    total_focal_associations += static_cast<double>(focal_indices.size());
    int n_target_hits = 0;
    for (int j = 0; j < n_targets; ++j) {
      if (target_hit[j]) ++n_target_hits;
    }
    total_focal_target_associations +=
      static_cast<double>(focal_indices.size()) * static_cast<double>(n_target_hits);

    for (std::size_t i = 0; i < focal_indices.size(); ++i) {
      const int origin_idx = focal_indices[i];
      tracks_counts[origin_idx] += 1;
      local_counts[origin_idx] += local_flag;
      for (int j = 0; j < n_targets; ++j) {
        if (target_hit[j]) {
          target_counts(origin_idx, j) += 1;
        }
      }
      if (!any_target && local_flag == 0) {
        leave_no_target_counts[origin_idx] += 1;
      }
    }
  }

  std::vector<int> penetration_target;
  std::vector<int> penetration_voxels;
  std::vector<double> penetration_counts;
  if (endpoint_analysis) {
    for (int target = 0; target < n_targets; ++target) {
      std::vector<std::pair<int, double> > bins(
        penetration_hist[target].begin(), penetration_hist[target].end()
      );
      std::sort(bins.begin(), bins.end(),
        [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
          return a.first < b.first;
        }
      );
      for (std::size_t i = 0; i < bins.size(); ++i) {
        penetration_target.push_back(target + 1);
        penetration_voxels.push_back(bins[i].first);
        penetration_counts.push_back(bins[i].second);
      }
    }
  }
  DataFrame penetration = DataFrame::create(
    _["target_index"] = IntegerVector(
      penetration_target.begin(), penetration_target.end()
    ),
    _["unique_target_voxels"] = IntegerVector(
      penetration_voxels.begin(), penetration_voxels.end()
    ),
    _["streamline_count"] = NumericVector(
      penetration_counts.begin(), penetration_counts.end()
    )
  );

  List endpoint = List::create(
    _["valid_streamlines"] = endpoint_valid_streamlines,
    _["endpoint_count"] = endpoint_count,
    _["endpoint_out_of_bounds"] = endpoint_out_of_bounds,
    _["endpoint_nonfinite"] = endpoint_nonfinite,
    _["endpoint_outside_one_voxel"] = endpoint_outside_one_voxel,
    _["endpoint_outside_more_than_one_voxel"] =
      endpoint_outside_more_than_one_voxel,
    _["endpoint_out_of_bounds_low"] = endpoint_out_of_bounds_low,
    _["endpoint_out_of_bounds_high"] = endpoint_out_of_bounds_high,
    _["endpoint_out_of_bounds_low_distance"] =
      endpoint_out_of_bounds_low_distance,
    _["endpoint_out_of_bounds_high_distance"] =
      endpoint_out_of_bounds_high_distance,
    _["grid_any_endpoint_inside"] = grid_any_endpoint_inside,
    _["grid_both_endpoints_inside"] = grid_both_endpoints_inside,
    _["context_endpoint_inside"] = context_endpoint_inside,
    _["context_any_endpoint_inside"] = context_any_endpoint_inside,
    _["context_both_endpoints_inside"] = context_both_endpoints_inside,
    _["focal_endpoint_streamlines_exact"] =
      focal_endpoint_streamlines_exact,
    _["focal_endpoint_streamlines_dilated1"] =
      focal_endpoint_streamlines_dilated1,
    _["paired_endpoint_streamlines_exact"] =
      paired_endpoint_streamlines_exact,
    _["paired_endpoint_streamlines_dilated1"] =
      paired_endpoint_streamlines_dilated1,
    _["paired_intersection_streamlines_exact"] =
      paired_intersection_streamlines_exact,
    _["paired_intersection_streamlines_dilated1"] =
      paired_intersection_streamlines_dilated1,
    _["focal_endpoint_target_intersection_exact"] =
      focal_endpoint_target_intersection_exact,
    _["focal_endpoint_target_intersection_dilated1"] =
      focal_endpoint_target_intersection_dilated1,
    _["class_exact"] = endpoint_class_exact,
    _["class_dilated1"] = endpoint_class_dilated1,
    _["penetration"] = penetration
  );
  SEXP endpoint_output = R_NilValue;
  if (endpoint_analysis) endpoint_output = endpoint;

  DataFrame local_search_assignments = DataFrame::create(
    _["chunk_streamline_index"] = IntegerVector(
      local_search_chunk_streamline_index.begin(),
      local_search_chunk_streamline_index.end()
    ),
    _["start_status_code"] = IntegerVector(
      local_search_start_label.begin(), local_search_start_label.end()
    ),
    _["end_status_code"] = IntegerVector(
      local_search_end_label.begin(), local_search_end_label.end()
    ),
    _["start_distance_mm"] = NumericVector(
      local_search_start_distance.begin(), local_search_start_distance.end()
    ),
    _["end_distance_mm"] = NumericVector(
      local_search_end_distance.begin(), local_search_end_distance.end()
    ),
    _["start_ambiguous"] = LogicalVector(
      local_search_start_ambiguous.begin(),
      local_search_start_ambiguous.end()
    ),
    _["end_ambiguous"] = LogicalVector(
      local_search_end_ambiguous.begin(), local_search_end_ambiguous.end()
    ),
    _["paired_target_index"] = IntegerVector(
      local_search_paired_target_index.begin(),
      local_search_paired_target_index.end()
    ),
    _["paired_focal_distance_mm"] = NumericVector(
      local_search_paired_focal_distance.begin(),
      local_search_paired_focal_distance.end()
    ),
    _["paired_target_distance_mm"] = NumericVector(
      local_search_paired_target_distance.begin(),
      local_search_paired_target_distance.end()
    ),
    _["paired_max_distance_mm"] = NumericVector(
      local_search_paired_max_distance.begin(),
      local_search_paired_max_distance.end()
    ),
    _["paired_target_intersects_full_path"] = LogicalVector(
      local_search_paired_full_path.begin(),
      local_search_paired_full_path.end()
    )
  );

  return List::create(
    _["tracks"] = tracks_counts,
    _["local"] = local_counts,
    _["leave_no_target"] = leave_no_target_counts,
    _["targets"] = target_counts,
    _["processed"] = processed,
    _["target_unique"] = target_unique_counts,
    _["full_path_target_unique"] = full_path_target_unique_counts,
	  _["local_search_assignments"] = local_search_assignments,
	    _["endpoint"] = endpoint_output,
    _["summary"] = NumericVector::create(
	      _["unique_focal_intersecting_streamlines"] = unique_focal_streamlines,
	      _["unique_local_streamlines"] = unique_local_streamlines,
	      _["unique_leave_no_target_streamlines"] = unique_leave_no_target_streamlines,
	      _["unique_focal_intersecting_with_target_streamlines"] = unique_with_target_streamlines,
	      _["unique_local_target_positive_streamlines"] = unique_local_target_positive_streamlines,
	      _["unique_local_target_negative_streamlines"] = unique_local_target_negative_streamlines,
	      _["unique_nonlocal_target_positive_streamlines"] = unique_nonlocal_target_positive_streamlines,
	      _["unique_nonlocal_target_negative_streamlines"] = unique_nonlocal_target_negative_streamlines,
	      _["unique_any_target_positive_streamlines"] = unique_any_target_positive_streamlines,
	      _["unique_any_target_negative_streamlines"] = unique_any_target_negative_streamlines,
	      _["total_focal_voxel_streamline_associations"] = total_focal_associations,
      _["total_focal_voxel_streamline_target_associations"] = total_focal_target_associations
    )
  );
}

// [[Rcpp::export]]
List ftq_qc_sample_chunk_cpp(List tracks,
                             NumericMatrix raw_to_voxel,
                             CharacterVector methods,
                             LogicalVector axis_flip,
                             IntegerVector vol_dim,
                             LogicalVector lut_f,
                             LogicalMatrix target_lut,
                             IntegerVector sample_quota,
                             IntegerVector current_counts,
                             std::string streamline_voxelization_mode,
                             std::string target_intersection_mode,
                             double local_search_radius_mm,
                             NumericMatrix local_search_voxel_to_world,
                             NumericVector local_search_voxel_radius) {
  if (raw_to_voxel.nrow() < 3 || raw_to_voxel.ncol() < 4) {
    stop("raw_to_voxel must have at least 3 rows and 4 columns.");
  }
  if (axis_flip.size() != 3 || vol_dim.size() < 3) {
    stop("axis_flip and vol_dim must have length 3.");
  }

  const std::vector<int> m_codes = method_codes(methods);
  const bool interpolate_segments = streamline_voxelization_mode == "segment_interpolation";
  if (streamline_voxelization_mode != "sample_points" && !interpolate_segments) {
    stop("Unsupported streamline_voxelization_mode in C++ QC sampler: " +
         streamline_voxelization_mode);
  }
  const TargetIntersectionMode target_mode = target_intersection_mode_code(
    target_intersection_mode, "C++ QC sampler"
  );
  if (target_mode == TARGET_LOCAL_SEARCH) {
    validate_local_search_parameters(
      local_search_radius_mm, local_search_voxel_to_world,
      local_search_voxel_radius, "C++ QC sampler"
    );
  }
  const int n_targets = target_lut.ncol();
  const int n_voxels = vol_dim[0] * vol_dim[1] * vol_dim[2];
  if (target_lut.nrow() != n_voxels || lut_f.size() != n_voxels) {
    stop("Focal and target lookup tables must match the reference grid.");
  }
  if (sample_quota.size() != n_targets || current_counts.size() != n_targets) {
    stop("sample_quota and current_counts must match target_lut columns.");
  }

  std::vector<int> counts(n_targets, 0);
  for (int j = 0; j < n_targets; ++j) counts[j] = current_counts[j];

  std::vector<int> sampled_track_indices;
  std::vector<int> sampled_target_indices;
  std::vector<int> target_seen_counts(n_targets, 0);
  int processed = 0;
  int focal_hits = 0;

  auto done = [&]() {
    for (int j = 0; j < n_targets; ++j) {
      if (sample_quota[j] > 0 && counts[j] < sample_quota[j]) return false;
    }
    return true;
  };

  for (R_xlen_t t = 0; t < tracks.size(); ++t) {
    if (done()) break;
    NumericMatrix coords = track_coords(tracks[t]);
    ++processed;
    if (coords.ncol() < 3 || coords.nrow() == 0) continue;

    std::vector<int> chosen_keys;
    bool focal_hit = false;

    for (std::size_t mi = 0; mi < m_codes.size(); ++mi) {
      std::vector<int> keys = voxelized_track_keys(
        coords, raw_to_voxel, m_codes[mi], axis_flip, vol_dim,
        interpolate_segments
      );
      bool method_focal_hit = false;

      for (std::size_t i = 0; i < keys.size(); ++i) {
        const int key = keys[i];
        if (key >= 0 && lut_f[key] == TRUE) method_focal_hit = true;
      }

      if (method_focal_hit) {
        focal_hit = true;
        chosen_keys.swap(keys);
        break;
      }
    }

    if (!focal_hit) continue;
    ++focal_hits;

    std::vector<int> full_path_target_hit(n_targets, 0);
    std::unordered_set<int> seen;
    seen.reserve(chosen_keys.size());
    for (std::size_t i = 0; i < chosen_keys.size(); ++i) {
      const int key = chosen_keys[i];
      if (key < 0) continue;
      if (!seen.insert(key).second) continue;
      for (int j = 0; j < n_targets; ++j) {
        if (target_lut(key, j) == TRUE) full_path_target_hit[j] = 1;
      }
    }
    std::vector<int> target_hit = primary_target_hits(
      target_mode, full_path_target_hit, chosen_keys, coords,
      raw_to_voxel, axis_flip, vol_dim, lut_f, target_lut,
      local_search_radius_mm, local_search_voxel_to_world,
      local_search_voxel_radius
    );

    int chosen_target = -1;
    int chosen_count = std::numeric_limits<int>::max();
    for (int j = 0; j < n_targets; ++j) {
      if (!target_hit[j]) continue;
      target_seen_counts[j] += 1;
      if (sample_quota[j] <= 0 || counts[j] >= sample_quota[j]) continue;
      if (counts[j] < chosen_count) {
        chosen_count = counts[j];
        chosen_target = j;
      }
    }
    if (chosen_target < 0) continue;

    sampled_track_indices.push_back(static_cast<int>(t) + 1);
    sampled_target_indices.push_back(chosen_target + 1);
    counts[chosen_target] += 1;
  }

  IntegerVector out_counts(n_targets);
  IntegerVector out_seen(n_targets);
  for (int j = 0; j < n_targets; ++j) {
    out_counts[j] = counts[j];
    out_seen[j] = target_seen_counts[j];
  }

  return List::create(
    _["track_indices"] = IntegerVector(sampled_track_indices.begin(), sampled_track_indices.end()),
    _["target_indices"] = IntegerVector(sampled_target_indices.begin(), sampled_target_indices.end()),
    _["counts"] = out_counts,
    _["target_seen"] = out_seen,
    _["processed"] = processed,
    _["focal_hits"] = focal_hits,
    _["complete"] = done()
  );
}
