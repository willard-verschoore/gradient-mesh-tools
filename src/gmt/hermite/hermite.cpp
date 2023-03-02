#include "gmt/hermite/hermite.hpp"

#include <algorithm>

namespace gmt::hermite
{

static CurveMatrix vec(float t)
{
  CurveMatrix row{{1.0f, t, t * t, t * t * t}};
  return row * PatchMatrix::N;
}

static CurveMatrix v_prime(float t)
{
  CurveMatrix row{{0.0f, 1.0f, 2.0f * t, 3.0f * t * t}};
  return row * PatchMatrix::N;
}

static CurveMatrix v_second(float t)
{
  CurveMatrix row{{0.0f, 0.0f, 2.0f, 6.0f * t}};
  return row * PatchMatrix::N;
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

Interpolant interpolate(const PatchMatrix& patch, float u, float v)
{
  return dot(vec(u) * patch, vec(v));
}

Interpolant parallel_derivative(const PatchMatrix& patch, float u, float v)
{
  return dot(vec(u) * patch, v_prime(v));
}

Interpolant orthogonal_derivative(const PatchMatrix& patch, float u, float v)
{
  return dot(v_prime(u) * patch, vec(v));
}

Interpolant mixed_derivative(const PatchMatrix& patch, float u, float v)
{
  return dot(v_prime(u) * patch, v_prime(v));
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

} // namespace gmt::hermite
