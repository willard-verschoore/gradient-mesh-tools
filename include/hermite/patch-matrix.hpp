#pragma once

#include <array>

#include "interpolant.hpp"

namespace hermite
{

struct PatchMatrix
{
  Interpolant data[4 * 4];

  // TODO: Consider making these local to hermite.cpp.
  static const PatchMatrix N;
  static const PatchMatrix N_INV;

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

} // namespace hermite
