#pragma once

#include <array>

#include "interpolant.hpp"

namespace gmt::hermite
{

struct PatchMatrix
{
  Interpolant data[4 * 4];

  // A Hermite patch with PatchMatrix G is given by:
  // H(u, v) = (1, u, u^2, u^3) * N * G * N^T * (1, v, v^2, v^3)^T.
  static const PatchMatrix N;
  static const PatchMatrix N_INV;

  // A Bezier patch with PatchMatrix P is given by:
  // B(u, v) = (1, u, u^2, u^3) * M * P * M^T * (1, v, v^2, v^3)^T.
  static const PatchMatrix M;
  static const PatchMatrix M_INV;

  constexpr PatchMatrix() = default;

  constexpr PatchMatrix(Interpolant const *values)
      : data{values[0],  values[1],  values[2],  values[3],
             values[4],  values[5],  values[6],  values[7],
             values[8],  values[9],  values[10], values[11],
             values[12], values[13], values[14], values[15]}
  {
  }

  constexpr PatchMatrix(std::array<Interpolant, 16> const &values)
      : data{values[0],  values[1],  values[2],  values[3],
             values[4],  values[5],  values[6],  values[7],
             values[8],  values[9],  values[10], values[11],
             values[12], values[13], values[14], values[15]}
  {
  }

  PatchMatrix &add(PatchMatrix const &other);
  PatchMatrix &add(Interpolant const &value);

  PatchMatrix &operator+=(PatchMatrix const &other);
  PatchMatrix &operator+=(Interpolant const &value);

  PatchMatrix &subtract(PatchMatrix const &other);
  PatchMatrix &subtract(Interpolant const &value);

  PatchMatrix &operator-=(PatchMatrix const &other);
  PatchMatrix &operator-=(Interpolant const &value);

  PatchMatrix &multiply(PatchMatrix const &other);
  PatchMatrix &multiply(Interpolant const &value);

  PatchMatrix &operator*=(PatchMatrix const &other);
  PatchMatrix &operator*=(Interpolant const &value);

  Interpolant &operator()(int row, int column);
  Interpolant const &operator()(int row, int column) const;

  PatchMatrix transposed() const;
  void transpose();

  PatchMatrix hermite() const;
  PatchMatrix bezier() const;
  void to_hermite();
  void to_bezier();
};

PatchMatrix add(PatchMatrix const &left, PatchMatrix const &right);
PatchMatrix add(PatchMatrix const &left, Interpolant const &right);
PatchMatrix add(Interpolant const &left, PatchMatrix const &right);

PatchMatrix operator+(PatchMatrix const &left, PatchMatrix const &right);
PatchMatrix operator+(PatchMatrix const &left, Interpolant const &right);
PatchMatrix operator+(Interpolant const &left, PatchMatrix const &right);

PatchMatrix subtract(PatchMatrix const &left, PatchMatrix const &right);
PatchMatrix subtract(PatchMatrix const &left, Interpolant const &right);
PatchMatrix subtract(Interpolant const &left, PatchMatrix const &right);

PatchMatrix operator-(PatchMatrix const &left, PatchMatrix const &right);
PatchMatrix operator-(PatchMatrix const &left, Interpolant const &right);
PatchMatrix operator-(Interpolant const &left, PatchMatrix const &right);

PatchMatrix multiply(PatchMatrix const &left, PatchMatrix const &right);
PatchMatrix multiply(PatchMatrix const &left, Interpolant const &right);
PatchMatrix multiply(Interpolant const &left, PatchMatrix const &right);

PatchMatrix operator*(PatchMatrix const &left, PatchMatrix const &right);
PatchMatrix operator*(PatchMatrix const &left, Interpolant const &right);
PatchMatrix operator*(Interpolant const &left, PatchMatrix const &right);

PatchMatrix negate(PatchMatrix const &patch);
PatchMatrix operator-(PatchMatrix const &patch);

} // namespace gmt::hermite
