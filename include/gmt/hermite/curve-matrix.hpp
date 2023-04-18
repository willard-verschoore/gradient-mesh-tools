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

  CurveMatrix hermite() const;
  CurveMatrix bezier() const;
  void to_hermite();
  void to_bezier();
};

Interpolant dot(CurveMatrix const &left, CurveMatrix const &right);

/// Checks whether \c left and \c right are equal.
/**
 * Two CurveMatrix objects are equal if each corresponding pair of elements is
 * equal.
 *
 * @param left,right The CurveMatrix objects to compare.
 * @return Whether \c left and \c right are equal.
 */
bool operator==(CurveMatrix const &left, CurveMatrix const &right);

/// Checks whether \c left and \c right are not equal.
/**
 * Two CurveMatrix objects are equal if each corresponding pair of elements is
 * equal. They are therefore unequal if there is at least one unequal pair.
 *
 * @param left,right The CurveMatrix objects to compare.
 * @return Whether \c left and \c right are not equal.
 */
bool operator!=(CurveMatrix const &left, CurveMatrix const &right);

} // namespace gmt::hermite
