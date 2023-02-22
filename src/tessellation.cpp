#include "tessellation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

using namespace hermite;

static float sample_distance(float max_deriv, float tolerance,
                             float view_length)
{
  auto displacement = 2.0f * max_deriv / view_length;
  return std::min(1.0f, std::sqrt(8.0f * tolerance / displacement));
}

static void compute_tess_levels(EdgeBoundary &boundary, float window_width,
                                float window_height, float viewport_width,
                                float viewport_height, float tolerance,
                                float max_tess, bool account_for_color)
{
  auto tx = tolerance / window_width;
  auto ty = tolerance / window_height;

  // Since color across two sides of the edge doesn't have to be contiguous,
  // color is not taken into account when computing the tessellation level,
  // but is instead given as a minimum tessellation level.
  auto minimum =
      account_for_color ? std::floor(std::pow(1.0f / tolerance, 2.0f)) : 1.0f;
  float maximum = std::floor(
      max_tess / std::max(boundary.segments.size(), boundary.num_siblings));
  for (auto &segment : boundary.segments)
  {
    auto max_deriv =
        max_second_parallel_derivative(boundary.curve, segment.interval);

    auto distance =
        std::min(sample_distance(max_deriv[0], tx, viewport_width),
                 sample_distance(max_deriv[1], ty, viewport_height));

    segment.tess_level =
        std::min(std::max(std::nearbyint(1.0f / distance), minimum), maximum);
  }
}

void cull_non_visible(std::vector<PatchRenderData> &patches,
                      float viewport_width, float viewport_height)
{
  // TODO: Reimplement.
  /*
  patches.erase(
      std::remove_if(patches.begin(), patches.end(),
                     [&](PatchRenderData const &patch) {
                       return !bounding_box(patch.matrix).intersects(viewport);
                     }),
      patches.end());
  */
}

void tessellate_patches(std::vector<PatchRenderData> &patches,
                        float window_width, float window_height,
                        float viewport_width, float viewport_height,
                        float tolerance, float max_hardware_tessellation,
                        bool account_for_color)
{
  cull_non_visible(patches, viewport_width, viewport_height);
  for (auto &p : patches)
  {
    compute_tess_levels(p.boundaries.top, window_width, window_height,
                        viewport_width, viewport_height, tolerance,
                        max_hardware_tessellation, account_for_color);
    compute_tess_levels(p.boundaries.left, window_width, window_height,
                        viewport_width, viewport_height, tolerance,
                        max_hardware_tessellation, account_for_color);
    compute_tess_levels(p.boundaries.bottom, window_width, window_height,
                        viewport_width, viewport_height, tolerance,
                        max_hardware_tessellation, account_for_color);
    compute_tess_levels(p.boundaries.right, window_width, window_height,
                        viewport_width, viewport_height, tolerance,
                        max_hardware_tessellation, account_for_color);
  }
}
