#include "draggable.hpp"

bool is_hovering(const QVector2D& cursor, const QVector2D& point,
                 const QRectF& viewport)
{
  return cursor.distanceToPoint(point) <= POINT_SIZE * viewport.height();
}
