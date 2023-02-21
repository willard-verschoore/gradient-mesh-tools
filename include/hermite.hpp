#pragma once

/**
 *  This file contains utility functions and data structures
 *  related to cubic Hermite splines and (bicubic) patches.
 */

#include <QGenericMatrix>
#include <QRectF>
#include <array>
#include <cstdint>

#include "interpolant.hpp"
#include "interval.hpp"

namespace hermite
{

/// Helper function to create a generic matrix from a std::array
template <std::size_t N, std::size_t M>
inline QGenericMatrix<N, M, Interpolant> matrix(
    const std::array<Interpolant, N * M>& arr)
{
  return QGenericMatrix<N, M, Interpolant>(arr.data());
}

using CurveMatrix = QGenericMatrix<4, 1, Interpolant>;
using PatchMatrix = QGenericMatrix<4, 4, Interpolant>;

/// Helper function to extract a specific row from a PatchMatrix.
inline CurveMatrix row(PatchMatrix mat, int row_index)
{
  return matrix<4, 1>({mat(row_index, 0), mat(row_index, 1), mat(row_index, 2),
                       mat(row_index, 3)});
}

static const PatchMatrix n = matrix<4, 4>({1.0f, 0.0f, 0.0f, 0.0f,    //
                                           0.0f, 1.0f, 0.0f, 0.0f,    //
                                           -3.0f, -2.0f, -1.0f, 3.0f, //
                                           2.0f, 1.0f, 1.0f, -2.0f});

static const PatchMatrix n_inv = matrix<4, 4>({1.0f, 0.0f, 0.0f, 0.0f, //
                                               0.0f, 1.0f, 0.0f, 0.0f, //
                                               0.0f, 1.0f, 2.0f, 3.0f, //
                                               1.0f, 1.0f, 1.0f, 1.0f});

/// Evaluate the given hermite curve at point `t`.
Interpolant interpolate(const CurveMatrix& mat, float t);
/// Evaluate the given hermite curve's first derivative at point `t`.
Interpolant parallel_derivative(const CurveMatrix& mat, float t);
/// Evaluate the given hermite curve's second derivative at point `t`.
Interpolant second_parallel_derivative(const CurveMatrix& mat, float t);
/**
 * Returns a vector containing the maximum absolute value
 * of the curve's second derivative for each component.
 * @param mat: the matrix describing the given curve.
 * @param interval: the interval in [0, 1] over which the maximum interval
 *      should be computed. [0, 1] is  the default.
 */
Interpolant max_second_parallel_derivative(const CurveMatrix& mat,
                                           const Interval& interval = {0.0f,
                                                                       1.0f});

/// Evaluate the given hermite patch at point `(u, v)`.
Interpolant interpolate(const PatchMatrix& mat, float u, float v);
/// Evaluate the given hermite patch's partial v-derivative at point `(u, v)`.
Interpolant parallel_derivative(const PatchMatrix& mat, float u, float v);
/// Evaluate the given hermite patch's partial u-derivative at point `(u, v)`.
Interpolant orthogonal_derivative(const PatchMatrix& mat, float u, float v);
/// Evaluate the given hermite patch's mixed uv-derivative at point `(u, v)`.
Interpolant mixed_derivative(const PatchMatrix& mat, float u, float v);

/// Returns the patch's AABB as a rectangle in world space.
QRectF bounding_box(const PatchMatrix& patch);
/// Returns an approximation of the point's distance to the given hermite curve.
float distance_to_curve(const Vector2& cursor, const CurveMatrix& curve);
} // namespace hermite
