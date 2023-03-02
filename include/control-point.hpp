#pragma once

#include "hermite/vector2.hpp"
#include "storage.hpp"

struct HalfEdge;

/// A user-editable control point. Each patch is linked to four control points.
struct ControlPoint
{
  /// Initialize the control point at the given coordinates.
  /**
   * @param coords: the control point's position in world space.
   */
  ControlPoint(hermite::Vector2 coords) : coords(coords) {}

  /// The control point's position in world space.
  hermite::Vector2 coords;
  /// Reference to the half-edge originating from this control point.
  Id<HalfEdge> edge;
};
