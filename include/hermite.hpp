#pragma once

/**
 *  This file contains utility functions and data structures
 *  related to cubic Hermite splines and (bicubic) patches.
 */

#include <array>
#include <cstdint>

#include "curve-matrix.hpp"
#include "interpolant.hpp"
#include "interval.hpp"
#include "patch-matrix.hpp"

namespace hermite
{

/// Evaluate the given hermite curve at point `t`.
Interpolant interpolate(const CurveMatrix& curve, float t);
/// Evaluate the given hermite curve's first derivative at point `t`.
Interpolant parallel_derivative(const CurveMatrix& curve, float t);
/// Evaluate the given hermite curve's second derivative at point `t`.
Interpolant second_parallel_derivative(const CurveMatrix& curve, float t);
/**
 * Returns a vector containing the maximum absolute value
 * of the curve's second derivative for each component.
 * @param curve: the matrix describing the given curve.
 * @param interval: the interval in [0, 1] over which the maximum interval
 *      should be computed. [0, 1] is  the default.
 */
Interpolant max_second_parallel_derivative(const CurveMatrix& curve,
                                           const Interval& interval = {0.0f,
                                                                       1.0f});

/// Evaluate the given hermite patch at point `(u, v)`.
Interpolant interpolate(const PatchMatrix& patch, float u, float v);
/// Evaluate the given hermite patch's partial v-derivative at point `(u, v)`.
Interpolant parallel_derivative(const PatchMatrix& patch, float u, float v);
/// Evaluate the given hermite patch's partial u-derivative at point `(u, v)`.
Interpolant orthogonal_derivative(const PatchMatrix& patch, float u, float v);
/// Evaluate the given hermite patch's mixed uv-derivative at point `(u, v)`.
Interpolant mixed_derivative(const PatchMatrix& patch, float u, float v);

/// Returns an approximation of the point's distance to the given hermite curve.
float distance_to_curve(const Vector2& cursor, const CurveMatrix& curve);

} // namespace hermite
