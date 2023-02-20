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

static void compute_tess_levels(EdgeBoundary& boundary,
                                const QSizeF& window_size,
                                const QRectF& viewport, float tolerance,
                                float max_tess, bool account_for_color)
{
  auto tx = tolerance / window_size.width();
  auto ty = tolerance / window_size.height();

  // Since color across two sides of the edge doesn't have to be contiguous,
  // color is not taken into account when computing the tessellation level,
  // but is instead given as a minimum tessellation level.
  auto minimum =
      account_for_color ? std::floor(std::pow(1.0f / tolerance, 2.0f)) : 1.0f;
  float maximum = std::floor(
      max_tess / std::max(boundary.segments.size(), boundary.num_siblings));
  for (auto& segment : boundary.segments)
  {
    auto max_deriv =
        max_second_parallel_derivative(boundary.curve, segment.interval);

    auto distance =
        std::min(sample_distance(max_deriv[0], tx, viewport.width()),
                 sample_distance(max_deriv[1], ty, viewport.height()));

    segment.tess_level =
        std::min(std::max(std::nearbyint(1.0f / distance), minimum), maximum);
  }
}

void cull_non_visible(std::vector<PatchRenderData>& patches,
                      const QRectF& viewport)
{
  patches.erase(
      std::remove_if(patches.begin(), patches.end(),
                     [&](const PatchRenderData& patch) {
                       return !bounding_box(patch.matrix).intersects(viewport);
                     }),
      patches.end());
}

void tessellate_patches(std::vector<PatchRenderData>& patches,
                        const QSizeF& window_size, const QRectF& viewport,
                        float tolerance, float max_hardware_tessellation,
                        bool account_for_color)
{
  cull_non_visible(patches, viewport);
  for (auto& p : patches)
  {
    compute_tess_levels(p.boundaries.top, window_size, viewport, tolerance,
                        max_hardware_tessellation, account_for_color);
    compute_tess_levels(p.boundaries.left, window_size, viewport, tolerance,
                        max_hardware_tessellation, account_for_color);
    compute_tess_levels(p.boundaries.bottom, window_size, viewport, tolerance,
                        max_hardware_tessellation, account_for_color);
    compute_tess_levels(p.boundaries.right, window_size, viewport, tolerance,
                        max_hardware_tessellation, account_for_color);
  }
}
