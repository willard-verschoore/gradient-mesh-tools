#include "hermite.hpp"

#include <algorithm>

namespace hermite
{

static CurveMatrix vec(float t)
{
  CurveMatrix row{{1.0f, t, t * t, t * t * t}};
  return dot(row, PatchMatrix::N);
}

static CurveMatrix v_prime(float t)
{
  CurveMatrix row{{0.0f, 1.0f, 2.0f * t, 3.0f * t * t}};
  return dot(row, PatchMatrix::N);
}

static CurveMatrix v_second(float t)
{
  CurveMatrix row{{0.0f, 0.0f, 2.0f, 6.0f * t}};
  return dot(row, PatchMatrix::N);
}

Interpolant interpolate(const CurveMatrix& curve, float t)
{
  return dot(curve, vec(t));
}

Interpolant parallel_derivative(const CurveMatrix& curve, float t)
{
  return dot(curve, v_prime(t));
}

Interpolant second_parallel_derivative(const CurveMatrix& curve, float t)
{
  return dot(curve, v_second(t));
}

Interpolant max_second_parallel_derivative(const CurveMatrix& curve,
                                           const Interval& interval)
{
  // Since second derivative is linear, the maximum and minimum values will be
  // at the endpoints.
  auto min = second_parallel_derivative(curve, interval.start);
  auto max = second_parallel_derivative(curve, interval.end);
  for (int i = 0; i < Interpolant::COMPONENTS; ++i)
    max[i] = std::max(std::abs(min[i]), std::abs(max[i]));
  return max * length(interval);
}

Interpolant interpolate(const PatchMatrix& patch, float u, float v)
{
  return dot(dot(vec(u), patch), vec(v));
}

Interpolant parallel_derivative(const PatchMatrix& patch, float u, float v)
{
  return dot(dot(vec(u), patch), v_prime(v));
}

Interpolant orthogonal_derivative(const PatchMatrix& patch, float u, float v)
{
  return dot(dot(v_prime(u), patch), vec(v));
}

Interpolant mixed_derivative(const PatchMatrix& curve, float u, float v)
{
  return dot(dot(v_prime(u), curve), v_prime(v));
}

static QRectF corner_bounds(const PatchMatrix& patch)
{
  auto top_left = QPointF(patch(0, 0).coords.x, patch(0, 0).coords.y);
  auto top_right = QPointF(patch(0, 3).coords.x, patch(0, 3).coords.y);
  auto bottom_left = QPointF(patch(3, 0).coords.x, patch(3, 0).coords.y);
  auto bottom_right = QPointF(patch(3, 3).coords.x, patch(3, 3).coords.y);

  return QRectF(top_left, bottom_right).united(QRectF(top_right, bottom_left));
}

static QPointF handle_position(const Interpolant& hermite_deriv,
                               const Interpolant& origin)
{
  Vector2 position = hermite_deriv.coords / 3.0f + origin.coords;
  return QPointF(position.x, position.y);
}

static QRectF handle_bounds(const PatchMatrix& patch)
{
  auto top_left = handle_position(patch(0, 1), patch(0, 0));
  auto top_right = handle_position(-patch(0, 2), patch(0, 3));
  auto bottom_left = handle_position(patch(3, 1), patch(3, 0));
  auto bottom_right = handle_position(-patch(3, 2), patch(3, 3));

  return QRectF(top_left, bottom_right).united(QRectF(top_right, bottom_left));
}

QRectF bounding_box(const PatchMatrix& patch)
{
  return corner_bounds(patch)
      .united(handle_bounds(patch))
      .united(handle_bounds(patch.transposed()));
}

float distance_to_curve(const Vector2& cursor, const CurveMatrix& curve)
{
  static constexpr std::size_t RESOLUTION = 5;
  auto dist = std::numeric_limits<float>::infinity();
  for (std::size_t step = 0; step < RESOLUTION; ++step)
  {
    auto point = interpolate(curve, float(step) / RESOLUTION).coords;
    dist = std::min(dist, distance(cursor, point));
  }
  return dist;
}

} // namespace hermite
