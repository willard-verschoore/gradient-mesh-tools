#include "gmt/hermite/patch-matrix.hpp"

#include <cstring>

namespace gmt::hermite
{

const PatchMatrix PatchMatrix::N{{1.0f, 0.0f, 0.0f, 0.0f,    //
                                  0.0f, 1.0f, 0.0f, 0.0f,    //
                                  -3.0f, -2.0f, -1.0f, 3.0f, //
                                  2.0f, 1.0f, 1.0f, -2.0f}};

const PatchMatrix PatchMatrix::N_INV{{1.0f, 0.0f, 0.0f, 0.0f, //
                                      0.0f, 1.0f, 0.0f, 0.0f, //
                                      0.0f, 1.0f, 2.0f, 3.0f, //
                                      1.0f, 1.0f, 1.0f, 1.0f}};

const PatchMatrix PatchMatrix::M{{1.0f, 0.0f, 0.0f, 0.0f,  //
                                  -3.0f, 3.0f, 0.0f, 0.0f, //
                                  3.0f, -6.0f, 3.0f, 0.0f, //
                                  -1.0f, 3.0f, -3.0f, 1.0f}};

const PatchMatrix PatchMatrix::M_INV{{1.0f, 0.0f, 0.0f, 0.0f,               //
                                      1.0f, 1.0f / 3.0f, 0.0f, 0.0f,        //
                                      1.0f, 2.0f / 3.0f, 1.0f / 3.0f, 0.0f, //
                                      1.0f, 1.0f, 1.0f, 1.0f}};

PatchMatrix &PatchMatrix::add(PatchMatrix const &other)
{
  for (int i = 0; i < 4 * 4; ++i) data[i] += other.data[i];
  return *this;
}

PatchMatrix &PatchMatrix::add(Interpolant const &value)
{
  for (int i = 0; i < 4 * 4; ++i) data[i] += value;
  return *this;
}

PatchMatrix &PatchMatrix::operator+=(PatchMatrix const &other)
{
  return add(other);
}

PatchMatrix &PatchMatrix::operator+=(Interpolant const &value)
{
  return add(value);
}

PatchMatrix &PatchMatrix::subtract(PatchMatrix const &other)
{
  for (int i = 0; i < 4 * 4; ++i) data[i] -= other.data[i];
  return *this;
}

PatchMatrix &PatchMatrix::subtract(Interpolant const &value)
{
  for (int i = 0; i < 4 * 4; ++i) data[i] -= value;
  return *this;
}

PatchMatrix &PatchMatrix::operator-=(PatchMatrix const &other)
{
  return subtract(other);
}

PatchMatrix &PatchMatrix::operator-=(Interpolant const &value)
{
  return subtract(value);
}

PatchMatrix &PatchMatrix::multiply(PatchMatrix const &other)
{
  Interpolant buffer[4 * 4];

  for (int i = 0; i < 4; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      buffer[j + i * 4] = data[0 + i * 4] * other.data[j + 0 * 4] +
                          data[1 + i * 4] * other.data[j + 1 * 4] +
                          data[2 + i * 4] * other.data[j + 2 * 4] +
                          data[3 + i * 4] * other.data[j + 3 * 4];
    }
  }

  memcpy(data, buffer, 4 * 4 * sizeof(Interpolant));
  return *this;
}

PatchMatrix &PatchMatrix::multiply(Interpolant const &value)
{
  for (int i = 0; i < 4 * 4; ++i) data[i] *= value;
  return *this;
}

PatchMatrix &PatchMatrix::operator*=(PatchMatrix const &other)
{
  return multiply(other);
}

PatchMatrix &PatchMatrix::operator*=(Interpolant const &value)
{
  return multiply(value);
}

Interpolant &PatchMatrix::operator()(int row, int column)
{
  return data[column + 4 * row];
}

Interpolant const &PatchMatrix::operator()(int row, int column) const
{
  return data[column + 4 * row];
}

PatchMatrix PatchMatrix::transposed() const
{
  return {{data[0], data[4], data[8], data[12],  //
           data[1], data[5], data[9], data[13],  //
           data[2], data[6], data[10], data[14], //
           data[3], data[7], data[11], data[15]}};
}

void PatchMatrix::transpose()
{
  Interpolant buffer[4 * 4] = {data[0], data[4], data[8],  data[12],
                               data[1], data[5], data[9],  data[13],
                               data[2], data[6], data[10], data[14],
                               data[3], data[7], data[11], data[15]};
  memcpy(data, buffer, 4 * 4 * sizeof(Interpolant));
}

PatchMatrix PatchMatrix::hermite() const
{
  return N_INV * M * (*this) * M.transposed() * N_INV.transposed();
}

PatchMatrix PatchMatrix::bezier() const
{
  return M_INV * N * (*this) * N.transposed() * M_INV.transposed();
}

void PatchMatrix::to_hermite()
{
  *this = N_INV * M * (*this) * M.transposed() * N_INV.transposed();
}

void PatchMatrix::to_bezier()
{
  *this = M_INV * N * (*this) * N.transposed() * M_INV.transposed();
}

PatchMatrix add(PatchMatrix const &left, PatchMatrix const &right)
{
  PatchMatrix copy{left};
  return copy.add(right);
}

PatchMatrix add(PatchMatrix const &left, Interpolant const &right)
{
  PatchMatrix copy{left};
  return copy.add(right);
}

PatchMatrix add(Interpolant const &left, PatchMatrix const &right)
{
  PatchMatrix copy{right};
  return copy.add(left);
}

PatchMatrix operator+(PatchMatrix const &left, PatchMatrix const &right)
{
  return add(left, right);
}

PatchMatrix operator+(PatchMatrix const &left, Interpolant const &right)
{
  return add(left, right);
}

PatchMatrix operator+(Interpolant const &left, PatchMatrix const &right)
{
  return add(left, right);
}

PatchMatrix subtract(PatchMatrix const &left, PatchMatrix const &right)
{
  PatchMatrix copy{left};
  return copy.subtract(right);
}

PatchMatrix subtract(PatchMatrix const &left, Interpolant const &right)
{
  PatchMatrix copy{left};
  return copy.subtract(right);
}

PatchMatrix subtract(Interpolant const &left, PatchMatrix const &right)
{
  PatchMatrix copy{right};
  return copy.subtract(left);
}

PatchMatrix operator-(PatchMatrix const &left, PatchMatrix const &right)
{
  return subtract(left, right);
}

PatchMatrix operator-(PatchMatrix const &left, Interpolant const &right)
{
  return subtract(left, right);
}

PatchMatrix operator-(Interpolant const &left, PatchMatrix const &right)
{
  return subtract(left, right);
}

PatchMatrix multiply(PatchMatrix const &left, PatchMatrix const &right)
{
  PatchMatrix copy{left};
  return copy.multiply(right);
}

PatchMatrix multiply(PatchMatrix const &left, Interpolant const &right)
{
  PatchMatrix copy{left};
  return copy.multiply(right);
}

PatchMatrix multiply(Interpolant const &left, PatchMatrix const &right)
{
  PatchMatrix copy{right};
  return copy.multiply(left);
}

PatchMatrix operator*(PatchMatrix const &left, PatchMatrix const &right)
{
  return multiply(left, right);
}

PatchMatrix operator*(PatchMatrix const &left, Interpolant const &right)
{
  return multiply(left, right);
}

PatchMatrix operator*(Interpolant const &left, PatchMatrix const &right)
{
  return multiply(left, right);
}

PatchMatrix negate(PatchMatrix const &patch)
{
  PatchMatrix copy{patch};
  for (int i = 0; i < 4 * 4; ++i) copy.data[i] = -patch.data[i];
  return copy;
}

PatchMatrix operator-(PatchMatrix const &patch) { return negate(patch); }

bool operator==(PatchMatrix const &left, PatchMatrix const &right)
{
  for (int i = 0; i < 4 * 4; ++i)
    if (left.data[i] != right.data[i]) return false;
  return true;
}

bool operator!=(PatchMatrix const &left, PatchMatrix const &right)
{
  return !(left == right);
}

std::ostream &operator<<(std::ostream &stream, PatchMatrix const &patch)
{
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 4; ++c)
      stream << '(' << r << ", " << c << ") " << patch(r, c) << '\n';
  return stream;
}
} // namespace gmt::hermite
