#pragma once

#include "storage.hpp"
#include "vector2.hpp"

struct HalfEdge;

/// A user-editable control point. Each patch is linked to four control points.
struct ControlPoint
{
  /// Initialize the control point at the given coordinates.
  /**
   * @param coords: the control point's position in world space.
   */
  ControlPoint(Vector2 coords) : coords(coords) {}

  /// The control point's position in world space.
  Vector2 coords;
  /// Reference to the half-edge originating from this control point.
  Id<HalfEdge> edge;
};
