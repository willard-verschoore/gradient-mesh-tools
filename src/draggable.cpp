#include "draggable.hpp"

bool is_hovering(const Vector2& cursor, const Vector2& point,
                 const QRectF& viewport)
{
  return distance(cursor, point) <= POINT_SIZE * viewport.height();
}
