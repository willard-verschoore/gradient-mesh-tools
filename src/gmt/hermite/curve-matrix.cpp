#include "gmt/hermite/curve-matrix.hpp"

namespace gmt::hermite
{

CurveMatrix::CurveMatrix(PatchMatrix const &patch, int row)
    : data{patch(row, 0), patch(row, 1), patch(row, 2), patch(row, 3)}
{
}

Interpolant &CurveMatrix::operator[](int index) { return data[index]; }

Interpolant const &CurveMatrix::operator[](int index) const
{
  return data[index];
}

Interpolant dot(CurveMatrix const &left, CurveMatrix const &right)
{
  return left[0] * right[0] + left[1] * right[1] + left[2] * right[2] +
         left[3] * right[3];
}

bool operator==(CurveMatrix const &left, CurveMatrix const &right)
{
  for (int i = 0; i < 4 * 4; ++i)
    if (left.data[i] != right.data[i]) return false;
  return true;
}

bool operator!=(CurveMatrix const &left, CurveMatrix const &right)
{
  return !(left == right);
}

} // namespace gmt::hermite
