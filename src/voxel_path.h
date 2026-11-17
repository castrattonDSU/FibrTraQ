#ifndef FIBRTRAQ_VOXEL_PATH_H
#define FIBRTRAQ_VOXEL_PATH_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
using namespace Rcpp;

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

} // namespace

#endif
