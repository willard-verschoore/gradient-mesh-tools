#pragma once

#include <QVector2D>

#include "draggable.hpp"
#include "storage.hpp"

struct HalfEdge;

/// A user-editable control point. Each patch is linked to four control points.
struct ControlPoint : public Draggable
{
  /// Initialize the control point at the given coordinates.
  /**
   * @param coords: the control point's position in world space.
   */
  ControlPoint(QVector2D coords) : coords(coords) {}

  void drag(const QVector2D& delta) override { coords += delta; }

  /// The control point's position in world space.
  QVector2D coords;
  /// Reference to the half-edge originating from this control point.
  Id<HalfEdge> edge;
};
