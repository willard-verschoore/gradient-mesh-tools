#include "gmt/hermite/hermite.hpp"

#include <catch2/catch.hpp>

using namespace gmt::hermite;

static CurveMatrix test_curve()
{
  CurveMatrix curve{{0.0f}};

  curve[0] = Interpolant{Vector2(0.0f, 0.0f), Vector3(1.0f, 0.0f, 0.0f)};
  curve[1] = Interpolant{Vector2(1.0f, 0.0f), Vector3(0.0f, 1.0f, 0.0f)};
  curve[2] = Interpolant{Vector2(0.0f, 1.0f), Vector3(0.0f, 0.0f, 1.0f)};
  curve[3] = Interpolant{Vector2(1.0f, 1.0f), Vector3(1.0f, 0.0f, 1.0f)};

  return curve;
}

static PatchMatrix test_patch()
{
  PatchMatrix patch{{0.0f}};

  patch(0, 0) = Interpolant{Vector2(0.0f, 0.0f), Vector3(1.0f, 0.0f, 0.0f)};
  patch(0, 3) = Interpolant{Vector2(1.0f, 0.0f), Vector3(0.0f, 1.0f, 0.0f)};
  patch(3, 0) = Interpolant{Vector2(0.0f, 1.0f), Vector3(0.0f, 0.0f, 1.0f)};
  patch(3, 3) = Interpolant{Vector2(1.0f, 1.0f), Vector3(1.0f, 0.0f, 1.0f)};

  patch(0, 1) = patch(0, 2) = Interpolant{Vector2(1.0f, 0.0f), Vector3()};
  patch(3, 1) = patch(3, 2) = Interpolant{Vector2(-1.0f, 0.0f), Vector3()};
  patch(1, 0) = patch(2, 0) = Interpolant{Vector2(0.0f, 1.0f), Vector3()};
  patch(1, 3) = patch(2, 3) = Interpolant{Vector2(0.0f, -1.0f), Vector3()};

  return patch;
}

static bool approx_eq(const Interpolant &a, const Interpolant &b)
{
  for (int i = 0; i < Interpolant::COMPONENTS; ++i)
  {
    // Always allow very small differences (useful if a = 0 or b = 0).
    if (std::abs(a[i] - b[i]) < 1e-4) continue;

    // Require smaller difference for small numbers and allow larger difference
    // for large numbers.
    if (std::abs(a[i] - b[i]) > std::max(std::abs(a[i]), std::abs(b[i])) * 1e-4)
      return false;
  }

  return true;
}

TEST_CASE("Surface interpolation at corner values", "[hermite]")
{
  const auto patch = test_patch();

  CHECK(approx_eq(patch(0, 0), interpolate(patch, 0.0f, 0.0f)));
  CHECK(approx_eq(patch(0, 3), interpolate(patch, 0.0f, 1.0f)));
  CHECK(approx_eq(patch(3, 0), interpolate(patch, 1.0f, 0.0f)));
  CHECK(approx_eq(patch(3, 3), interpolate(patch, 1.0f, 1.0f)));
}

TEST_CASE("Derivatives at corners are equal to handles", "[hermite]")
{
  const auto patch = test_patch();

  CHECK(approx_eq(patch(1, 0), orthogonal_derivative(patch, 0.0f, 0.0f)));
  CHECK(approx_eq(patch(2, 0), orthogonal_derivative(patch, 1.0f, 0.0f)));
  CHECK(approx_eq(patch(1, 3), orthogonal_derivative(patch, 0.0f, 1.0f)));
  CHECK(approx_eq(patch(2, 3), orthogonal_derivative(patch, 1.0f, 1.0f)));

  CHECK(approx_eq(patch(0, 1), parallel_derivative(patch, 0.0f, 0.0f)));
  CHECK(approx_eq(patch(0, 2), parallel_derivative(patch, 0.0f, 1.0f)));
  CHECK(approx_eq(patch(3, 1), parallel_derivative(patch, 1.0f, 0.0f)));
  CHECK(approx_eq(patch(3, 2), parallel_derivative(patch, 1.0f, 1.0f)));
}

TEST_CASE("Mixed derivatives at corners are equal to twist vectors",
          "[hermite]")
{
  const auto patch = test_patch();

  CHECK(approx_eq(patch(1, 1), mixed_derivative(patch, 0.0f, 0.0f)));
  CHECK(approx_eq(patch(2, 1), mixed_derivative(patch, 1.0f, 0.0f)));
  CHECK(approx_eq(patch(1, 2), mixed_derivative(patch, 0.0f, 1.0f)));
  CHECK(approx_eq(patch(2, 2), mixed_derivative(patch, 1.0f, 1.0f)));
}

TEST_CASE("Surface interpolation at non-corner values", "[hermite]")
{
  const auto patch = test_patch();

  auto top = Interpolant{Vector2(0.5f, 0.0f), Vector3(0.5f, 0.5f, 0.0f)};
  auto bottom = Interpolant{Vector2(0.5f, 1.0f), Vector3(0.5f, 0.0f, 1.0f)};
  auto left = Interpolant{Vector2(0.0f, 0.5f), Vector3(0.5f, 0.0f, 0.5f)};
  auto right = Interpolant{Vector2(1.0f, 0.5f), Vector3(0.5f, 0.5f, 0.5f)};
  auto center = Interpolant{Vector2(0.5f, 0.5f), Vector3(0.5f, 0.25f, 0.5f)};

  CHECK(approx_eq(top, interpolate(patch, 0.0f, 0.5f)));
  CHECK(approx_eq(bottom, interpolate(patch, 1.0f, 0.5f)));
  CHECK(approx_eq(left, interpolate(patch, 0.5f, 0.0f)));
  CHECK(approx_eq(right, interpolate(patch, 0.5f, 1.0f)));
  CHECK(approx_eq(center, interpolate(patch, 0.5f, 0.5f)));
}

TEST_CASE("Curve conversion to and from Bezier representation", "[hermite]")
{
  const auto curve = test_curve();

  auto in_place = test_curve();
  in_place.to_bezier();
  in_place.to_hermite();
  auto copy = curve.bezier().hermite();

  for (int i = 0; i < 4; ++i)
  {
    CHECK(approx_eq(curve[i], in_place[i]));
    CHECK(approx_eq(curve[i], copy[i]));
  }
}

TEST_CASE("Patch conversion to and from Bezier representation", "[hermite]")
{
  const auto patch = test_patch();

  auto in_place = test_patch();
  in_place.to_bezier();
  in_place.to_hermite();
  auto copy = patch.bezier().hermite();

  for (int i = 0; i < 16; ++i)
  {
    CHECK(approx_eq(patch.data[i], in_place.data[i]));
    CHECK(approx_eq(patch.data[i], copy.data[i]));
  }
}
