#pragma once

#include <QRectF>

#include "hermite.hpp"
#include "interval.hpp"
#include "patch.hpp"
#include "quad.hpp"

/// The data needed to render a patch boundary.
/**
 * and the locations of its abutting t-junctions
 * (useful for crack-free tessellation).
 */
struct EdgeBoundary
{
  struct Segment
  {
    Interval interval;
    float tess_level = 0.0f;
    /**
     * Unused field used as padding for the C++ struct definition
     * to align with the GLSL struct definition.
     */
    float _padding = 0.0f;
  };
  /// The boundary curve's control matrix in Hermite form.
  hermite::CurveMatrix curve;
  std::vector<Segment> segments;
  std::size_t num_siblings;
};

/// The data needed to render a patch and all four of its boundaries.
struct PatchRenderData
{
  PatchRenderData(hermite::PatchMatrix matrix, Quad<EdgeBoundary> boundaries)
      : matrix(std::move(matrix)), boundaries(std::move(boundaries))
  {
  }

  /// The boundary curve's control matrix in Hermite form.
  hermite::PatchMatrix matrix;
  Quad<EdgeBoundary> boundaries;
};

/**
 * Compute the tessellation level for all of the given edges,
 * updating all segments' tessellation levels.
 */
void tessellate_patches(std::vector<PatchRenderData>& patches,
                        const QSizeF& window_size, const QRectF& viewport,
                        float tolerance, float max_hardware_tessellation,
                        bool account_for_color);
