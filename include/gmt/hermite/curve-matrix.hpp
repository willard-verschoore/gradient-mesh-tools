#pragma once

#include "patch-matrix.hpp"

namespace gmt::hermite
{

struct CurveMatrix
{
  Interpolant data[4];

  constexpr CurveMatrix() = default;

  constexpr CurveMatrix(Interpolant const *values)
      : data{values[0], values[1], values[2], values[3]}
  {
  }

  constexpr CurveMatrix(std::array<Interpolant, 4> const &values)
      : data{values[0], values[1], values[2], values[3]}
  {
  }

  CurveMatrix(PatchMatrix const &patch, int row);

  Interpolant &operator[](int index);
  Interpolant const &operator[](int index) const;
};

Interpolant dot(CurveMatrix const &left, CurveMatrix const &right);

} // namespace gmt::hermite
