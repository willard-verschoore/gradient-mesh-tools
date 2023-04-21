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

/// Writes the given curve to the given stream.
/**
 * Will write the curve's four elements separately. Each element will be written
 * on it's own line with a number indicating which element is being printed.
 *
 * @param stream The output stream to write to.
 * @param vector The CurveMatrix object to write.
 * @return A reference to \c stream.
 */
std::ostream &operator<<(std::ostream &stream, CurveMatrix const &curve);

} // namespace gmt::hermite
