#pragma once

#include "hermite/hermite.hpp"
#include "interval.hpp"
#include "patch.hpp"

/// The data needed to render one edge of a patch boundary.
/**
 * This also specifies the locations of an edge's abutting t-junctions, which is
 * is useful information for crack-free tessellation.
 */
struct EdgeBoundary
{
  struct Segment
  {
    Interval interval;
    float tess_level = 0.0f;

    /**
     * Unused field used as padding for the C++ struct definition to align with
     * the GLSL struct definition.
     */
    float _padding = 0.0f;
  };

  /// The boundary curve's control matrix in Hermite form.
  hermite::CurveMatrix curve;
  std::vector<Segment> segments;
  std::size_t num_siblings;
};

/// The collection of four EdgeBoundary objects making up a patch boundary.
struct PatchBoundary
{
  EdgeBoundary top, left, bottom, right;
};

/// The data needed to render a patch and its entire boundary.
struct PatchRenderData
{
  PatchRenderData(hermite::PatchMatrix const &matrix,
                  PatchBoundary const &boundary)
      : matrix(matrix), boundary(boundary)
  {
  }

  /// The patch's matrix in Hermite form.
  hermite::PatchMatrix matrix;

  /// The patch boundary made up of four edge boundaries.
  PatchBoundary boundary;
};

/**
 * Compute the tessellation level for all of the given edges, updating all
 * segments' tessellation levels.
 */
void tessellate_patches(std::vector<PatchRenderData> &patches,
                        float window_width, float window_height,
                        float viewport_width, float viewport_height,
                        float tolerance, float max_hardware_tessellation,
                        bool account_for_color);
