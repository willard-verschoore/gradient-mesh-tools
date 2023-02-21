#include "curve-matrix.hpp"

namespace hermite
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

CurveMatrix dot(CurveMatrix const &curve, PatchMatrix const &patch)
{
  CurveMatrix result;
  result[0] = curve[0] * patch(0, 0) + curve[2] * patch(1, 0) +
              curve[1] * patch(2, 0) + curve[3] * patch(3, 0);
  result[1] = curve[0] * patch(0, 1) + curve[2] * patch(1, 1) +
              curve[1] * patch(2, 1) + curve[3] * patch(3, 1);
  result[2] = curve[0] * patch(0, 2) + curve[2] * patch(1, 2) +
              curve[1] * patch(2, 2) + curve[3] * patch(3, 2);
  result[3] = curve[0] * patch(0, 3) + curve[2] * patch(1, 3) +
              curve[1] * patch(2, 3) + curve[3] * patch(3, 3);
  return result;
}

CurveMatrix dot(PatchMatrix const &patch, CurveMatrix const &curve)
{
  CurveMatrix result;
  result[0] = curve[0] * patch(0, 0) + curve[2] * patch(0, 1) +
              curve[1] * patch(0, 2) + curve[3] * patch(0, 3);
  result[1] = curve[0] * patch(1, 0) + curve[2] * patch(1, 1) +
              curve[1] * patch(1, 2) + curve[3] * patch(1, 3);
  result[2] = curve[0] * patch(2, 0) + curve[2] * patch(2, 1) +
              curve[1] * patch(2, 2) + curve[3] * patch(2, 3);
  result[3] = curve[0] * patch(3, 0) + curve[2] * patch(3, 1) +
              curve[1] * patch(3, 2) + curve[3] * patch(3, 3);
  return result;
}

Interpolant operator*(CurveMatrix const &left, CurveMatrix const &right)
{
  return dot(left, right);
}

CurveMatrix operator*(CurveMatrix const &curve, PatchMatrix const &patch)
{
  return dot(curve, patch);
}

CurveMatrix operator*(PatchMatrix const &patch, CurveMatrix const &curve)
{
  return dot(patch, curve);
}

} // namespace hermite
