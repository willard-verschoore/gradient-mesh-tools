#pragma once

/**
 *  This file contains utility functions and data structures
 *  related to cubic Hermite splines and (bicubic) patches.
 */

#include <array>
#include <cstdint>

#include "curve-matrix.hpp"
#include "interpolant.hpp"
#include "patch-matrix.hpp"
#include "vector2.hpp"
#include "vector3.hpp"

namespace gmt::hermite
{

/// Returns the matrix multiplication of \c curve and \c patch.
/**
 * Performs the matrix multiplication \c curve * \c patch where \c curve is 1x4
 * and \c patch is 4x4. This means that \c curve is interpreted as a row vector.
 *
 * @param curve The CurveMatrix object to multiply.
 * @param patch The PatchMatrix object to multiply.
 * @return The CurveMatrix resulting from matrix multiplication of \c curve and
 * \c patch.
 */
CurveMatrix multiply(CurveMatrix const& curve, PatchMatrix const& patch);

/// Returns the matrix multiplication of \c patch and \c curve.
/**
 * Performs the matrix multiplication \c patch * \c curve where \c patch is 4x4
 * and \c curve is 4x1. This means that \c curve is interpreted as a column
 * vector.
 *
 * @param patch The PatchMatrix object to multiply.
 * @param curve The CurveMatrix object to multiply.
 * @return The CurveMatrix resulting from matrix multiplication of \c patch and
 * \c curve.
 */
CurveMatrix multiply(PatchMatrix const& patch, CurveMatrix const& curve);

/// Returns the matrix multiplication of \c curve and \c patch.
/**
 * Performs the matrix multiplication \c curve * \c patch where \c curve is 1x4
 * and \c patch is 4x4. This means that \c curve is interpreted as a row vector.
 *
 * @param curve The CurveMatrix object to multiply.
 * @param patch The PatchMatrix object to multiply.
 * @return The CurveMatrix resulting from matrix multiplication of \c curve and
 * \c patch.
 */
CurveMatrix operator*(CurveMatrix const& curve, PatchMatrix const& patch);

/// Returns the matrix multiplication of \c patch and \c curve.
/**
 * Performs the matrix multiplication \c patch * \c curve where \c patch is 4x4
 * and \c curve is 4x1. This means that \c curve is interpreted as a column
 * vector.
 *
 * @param patch The PatchMatrix object to multiply.
 * @param curve The CurveMatrix object to multiply.
 * @return The CurveMatrix resulting from matrix multiplication of \c patch and
 * \c curve.
 */
CurveMatrix operator*(PatchMatrix const& patch, CurveMatrix const& curve);

/// Evaluate the given hermite curve at point `t`.
Interpolant interpolate(const CurveMatrix& curve, float t);
/// Evaluate the given hermite curve's first derivative at point `t`.
Interpolant parallel_derivative(const CurveMatrix& curve, float t);
/// Evaluate the given hermite curve's second derivative at point `t`.
Interpolant second_parallel_derivative(const CurveMatrix& curve, float t);
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

} // namespace gmt::hermite
