#include "gmt/hermite/interpolant.hpp"

namespace gmt::hermite
{

float& Interpolant::operator[](int index)
{
  if (index < 2)
    return coords[index];
  else
    return color[index - 2];
}

float Interpolant::operator[](int index) const
{
  if (index < 2)
    return coords[index];
  else
    return color[index - 2];
}

Interpolant operator-(Interpolant i)
{
  i.coords = -i.coords;
  i.color = -i.color;
  return i;
}

Interpolant& operator/=(Interpolant& left, float right)
{
  left.coords /= right;
  left.color /= right;
  return left;
}

Interpolant& operator*=(Interpolant& left, float right)
{
  left.coords *= right;
  left.color *= right;
  return left;
}

Interpolant& operator*=(Interpolant& left, const Interpolant& right)
{
  left.coords *= right.coords;
  left.color *= right.color;
  return left;
}

Interpolant& operator+=(Interpolant& left, const Interpolant& right)
{
  left.coords += right.coords;
  left.color += right.color;
  return left;
}

Interpolant& operator-=(Interpolant& left, const Interpolant& right)
{
  left.coords -= right.coords;
  left.color -= right.color;
  return left;
}

bool is_zero(const Interpolant& i, float coords_epsilon, float color_epsilon)
{
  for (int k = 0; k < 2; ++k)
    if (std::abs(i.coords[k]) > coords_epsilon) return false;
  for (int k = 0; k < 3; ++k)
    if (std::abs(i.color[k]) > color_epsilon) return false;
  return true;
}

float length(const Interpolant& i)
{
  float size = 0;
  for (int k = 0; k < i.COMPONENTS; ++k) size += pow(i[k], 2);
  return sqrt(size);
}

} // namespace gmt::hermite
