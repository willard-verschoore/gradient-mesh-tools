#pragma once

#include <QRectF>
#include <QVector2D>

/// Abstract class for user-moveable controls (e.g. control points and handles)
/**
 * An abstract class for something which can be dragged in the user interface,
 * such as control points and tangent handles.
 */
class Draggable
{
 public:
  virtual ~Draggable() {}

  /// Move this object a distance of `delta` from its current position.
  /**
   * @param delta: the vector representing the difference between the object's
   * current and intended new position.
   */
  virtual void drag(const QVector2D& delta) = 0;
};

/**
 * The size of a control point as a fraction of the screen height,
 * to be used to determine whether the cursor is hovering a Draggable.
 */
inline float POINT_SIZE = 0.02f;

/// Determines if the cursor is hovering above the given point.
/**
 * @param cursor: the cursor position in world space.
 * @param point: the target point's position in world space.
 * @param viewport: the current viewport as a world-space rectangle.
 */
bool is_hovering(const QVector2D& cursor, const QVector2D& point,
                 const QRectF& viewport);
