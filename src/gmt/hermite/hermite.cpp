#include "gmt/hermite/hermite.hpp"

#include <algorithm>

namespace gmt::hermite
{

CurveMatrix multiply(CurveMatrix const &curve, PatchMatrix const &patch)
{
  CurveMatrix result;
  result[0] = curve[0] * patch(0, 0) + curve[1] * patch(1, 0) +
              curve[2] * patch(2, 0) + curve[3] * patch(3, 0);
  result[1] = curve[0] * patch(0, 1) + curve[1] * patch(1, 1) +
              curve[2] * patch(2, 1) + curve[3] * patch(3, 1);
  result[2] = curve[0] * patch(0, 2) + curve[1] * patch(1, 2) +
              curve[2] * patch(2, 2) + curve[3] * patch(3, 2);
  result[3] = curve[0] * patch(0, 3) + curve[1] * patch(1, 3) +
              curve[2] * patch(2, 3) + curve[3] * patch(3, 3);
  return result;
}

CurveMatrix multiply(PatchMatrix const &patch, CurveMatrix const &curve)
{
  CurveMatrix result;
  result[0] = curve[0] * patch(0, 0) + curve[1] * patch(0, 1) +
              curve[2] * patch(0, 2) + curve[3] * patch(0, 3);
  result[1] = curve[0] * patch(1, 0) + curve[1] * patch(1, 1) +
              curve[2] * patch(1, 2) + curve[3] * patch(1, 3);
  result[2] = curve[0] * patch(2, 0) + curve[1] * patch(2, 1) +
              curve[2] * patch(2, 2) + curve[3] * patch(2, 3);
  result[3] = curve[0] * patch(3, 0) + curve[1] * patch(3, 1) +
              curve[2] * patch(3, 2) + curve[3] * patch(3, 3);
  return result;
}

CurveMatrix operator*(CurveMatrix const &curve, PatchMatrix const &patch)
{
  return multiply(curve, patch);
}

CurveMatrix operator*(PatchMatrix const &patch, CurveMatrix const &curve)
{
  return multiply(patch, curve);
}

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

Interpolant interpolate(CurveMatrix const &curve, float t)
{
  return dot(curve, vec(t));
}

Interpolant parallel_derivative(CurveMatrix const &curve, float t)
{
  return dot(curve, v_prime(t));
}

Interpolant second_parallel_derivative(CurveMatrix const &curve, float t)
{
  return dot(curve, v_second(t));
}

Interpolant interpolate(PatchMatrix const &patch, float u, float v)
{
  return dot(vec(u) * patch, vec(v));
}

Interpolant parallel_derivative(PatchMatrix const &patch, float u, float v)
{
  return dot(vec(u) * patch, v_prime(v));
}

Interpolant orthogonal_derivative(PatchMatrix const &patch, float u, float v)
{
  return dot(v_prime(u) * patch, vec(v));
}

Interpolant mixed_derivative(PatchMatrix const &patch, float u, float v)
{
  return dot(v_prime(u) * patch, v_prime(v));
}

float distance_to_curve(Vector2 const &cursor, CurveMatrix const &curve)
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
