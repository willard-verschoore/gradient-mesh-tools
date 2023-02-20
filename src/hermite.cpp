#include "hermite.hpp"

#include <algorithm>

namespace hermite
{

template <std::size_t I>
static inline CurveMatrix extract_row(const PatchMatrix& mat)
{
  std::array<Interpolant, 8> arr;
  mat.copyDataTo(arr.data());
  return CurveMatrix(arr.data() + I * 4);
}

static inline Interpolant dot(const CurveMatrix& a, const CurveMatrix& b)
{
  return (a * b.transposed())(0, 0);
}

static CurveMatrix vec(float t)
{
  auto row = matrix<4, 1>({
      1.0f,
      t,
      t * t,
      t * t * t,
  });
  return row * n;
}

static CurveMatrix v_prime(float t)
{
  auto row = matrix<4, 1>({
      0.0f,
      1.0f,
      2.0f * t,
      3.0f * t * t,
  });
  return row * n;
}
static CurveMatrix v_second(float t)
{
  auto row = matrix<4, 1>({
      0.0f,
      0.0f,
      2.0f,
      6.0f * t,
  });
  return row * n;
}

Interpolant interpolate(const CurveMatrix& mat, float t)
{
  return dot(mat, vec(t));
}
Interpolant parallel_derivative(const CurveMatrix& mat, float t)
{
  return dot(mat, v_prime(t));
}
Interpolant second_parallel_derivative(const CurveMatrix& mat, float t)
{
  return dot(mat, v_second(t));
}

Interpolant max_second_parallel_derivative(const CurveMatrix& mat,
                                           const Interval& interval)
{
  // Since second derivative is linear, the maximum and minimum values will be
  // at the endpoints.
  auto min = second_parallel_derivative(mat, interval.start);
  auto max = second_parallel_derivative(mat, interval.end);
  for (int i = 0; i < Interpolant::COMPONENTS; ++i)
  {
    max[i] = std::max(std::abs(min[i]), std::abs(max[i]));
  }
  return max * length(interval);
}

Interpolant interpolate(const PatchMatrix& mat, float u, float v)
{
  return dot(vec(u) * mat, vec(v));
}
Interpolant parallel_derivative(const PatchMatrix& mat, float u, float v)
{
  return dot(vec(u) * mat, v_prime(v));
}
Interpolant orthogonal_derivative(const PatchMatrix& mat, float u, float v)
{
  return dot(v_prime(u) * mat, vec(v));
}
Interpolant mixed_derivative(const PatchMatrix& mat, float u, float v)
{
  return dot(v_prime(u) * mat, v_prime(v));
}

static QRectF corner_bounds(const PatchMatrix& patch)
{
  auto top_left = patch(0, 0).coords.toPointF();
  auto top_right = patch(0, 3).coords.toPointF();
  auto bottom_left = patch(3, 0).coords.toPointF();
  auto bottom_right = patch(3, 3).coords.toPointF();

  return QRectF(top_left, bottom_right).united(QRectF(top_right, bottom_left));
}

static QPointF handle_position(const Interpolant& hermite_deriv,
                               const Interpolant& origin)
{
  return (hermite_deriv.coords / 3.0f + origin.coords).toPointF();
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

float distance_to_curve(const QVector2D& cursor, const CurveMatrix& curve)
{
  static constexpr std::size_t RESOLUTION = 5;
  auto dist = std::numeric_limits<float>::infinity();
  for (std::size_t step = 0; step < RESOLUTION; ++step)
  {
    auto point = interpolate(curve, float(step) / RESOLUTION).coords;
    dist = std::min(dist, cursor.distanceToPoint(point));
  }
  return dist;
}
} // namespace hermite
