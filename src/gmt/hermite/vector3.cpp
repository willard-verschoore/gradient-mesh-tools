#include "gmt/hermite/vector3.hpp"

#include <cmath>

namespace gmt::hermite
{

Vector3 &Vector3::add(Vector3 const &other)
{
  x += other.x;
  y += other.y;
  z += other.z;

  return *this;
}

Vector3 &Vector3::add(float value)
{
  x += value;
  y += value;
  z += value;

  return *this;
}

Vector3 &Vector3::operator+=(Vector3 const &other) { return add(other); }

Vector3 &Vector3::operator+=(float value) { return add(value); }

Vector3 add(Vector3 const &left, Vector3 const &right)
{
  Vector3 copy{left};
  return copy.add(right);
}

Vector3 add(Vector3 const &left, float right)
{
  Vector3 copy{left};
  return copy.add(right);
}

Vector3 add(float left, Vector3 const &right)
{
  Vector3 copy{right};
  return copy.add(left);
}

Vector3 operator+(Vector3 const &left, Vector3 const &right)
{
  return add(left, right);
}

Vector3 operator+(Vector3 const &left, float right) { return add(left, right); }

Vector3 operator+(float left, Vector3 const &right) { return add(left, right); }

Vector3 &Vector3::subtract(Vector3 const &other)
{
  x -= other.x;
  y -= other.y;
  z -= other.z;

  return *this;
}

Vector3 &Vector3::subtract(float value)
{
  x -= value;
  y -= value;
  z -= value;

  return *this;
}

Vector3 &Vector3::operator-=(Vector3 const &other) { return subtract(other); }

Vector3 &Vector3::operator-=(float value) { return subtract(value); }

Vector3 subtract(Vector3 const &left, Vector3 const &right)
{
  Vector3 copy{left};
  return copy.subtract(right);
}

Vector3 subtract(Vector3 const &left, float right)
{
  Vector3 copy{left};
  return copy.subtract(right);
}

Vector3 operator-(Vector3 const &left, Vector3 const &right)
{
  return subtract(left, right);
}

Vector3 operator-(Vector3 const &left, float right)
{
  return subtract(left, right);
}

Vector3 &Vector3::multiply(Vector3 const &other)
{
  x *= other.x;
  y *= other.y;
  z *= other.z;

  return *this;
}

Vector3 &Vector3::multiply(float value)
{
  x *= value;
  y *= value;
  z *= value;

  return *this;
}

Vector3 &Vector3::operator*=(Vector3 const &other) { return multiply(other); }

Vector3 &Vector3::operator*=(float value) { return multiply(value); }

Vector3 multiply(Vector3 const &left, Vector3 const &right)
{
  Vector3 copy{left};
  return copy.multiply(right);
}

Vector3 multiply(Vector3 const &left, float right)
{
  Vector3 copy{left};
  return copy.multiply(right);
}

Vector3 multiply(float left, Vector3 const &right)
{
  Vector3 copy{right};
  return copy.multiply(left);
}

Vector3 operator*(Vector3 const &left, Vector3 const &right)
{
  return multiply(left, right);
}

Vector3 operator*(Vector3 const &left, float right)
{
  return multiply(left, right);
}

Vector3 operator*(float left, Vector3 const &right)
{
  return multiply(left, right);
}

Vector3 &Vector3::divide(Vector3 const &other)
{
  x /= other.x;
  y /= other.y;
  z /= other.z;

  return *this;
}

Vector3 &Vector3::divide(float value)
{
  float inverse = 1.0f / value;

  x *= inverse;
  y *= inverse;
  z *= inverse;

  return *this;
}

Vector3 &Vector3::operator/=(Vector3 const &other) { return divide(other); }

Vector3 &Vector3::operator/=(float value) { return divide(value); }

Vector3 divide(Vector3 const &left, Vector3 const &right)
{
  Vector3 copy{left};
  return copy.divide(right);
}

Vector3 divide(Vector3 const &left, float right)
{
  Vector3 copy{left};
  return copy.divide(right);
}

Vector3 operator/(Vector3 const &left, Vector3 const &right)
{
  return divide(left, right);
}

Vector3 operator/(Vector3 const &left, float right)
{
  return divide(left, right);
}

float &Vector3::operator[](int index) { return data[index]; }

float Vector3::operator[](int index) const { return data[index]; }

Vector3 negate(Vector3 const &value)
{
  Vector3 copy{value};

  copy.x = -copy.x;
  copy.y = -copy.y;
  copy.z = -copy.z;

  return copy;
}

Vector3 operator-(Vector3 const &value) { return negate(value); }

float Vector3::length() const { return sqrtf(length_squared()); }

float Vector3::length_squared() const { return x * x + y * y + z * z; }

Vector3 Vector3::normalized() const { return (*this) / length(); }

void Vector3::normalize()
{
  float inverse = 1.0f / length();

  x *= inverse;
  y *= inverse;
  z *= inverse;
}

float dot(Vector3 const &left, Vector3 const &right)
{
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

Vector3 cross(Vector3 const &left, Vector3 const &right)
{
  return Vector3(left.y * right.z - left.z * right.y,
                 left.z * right.x - left.x * right.z,
                 left.x * right.y - left.y * right.x);
}

float distance(Vector3 const &left, Vector3 const &right)
{
  return (left - right).length();
}

float distance_squared(Vector3 const &left, Vector3 const &right)
{
  return (left - right).length_squared();
}

std::ostream &operator<<(std::ostream &stream, Vector3 const &vector)
{
  return stream << "(" << vector.x << ", " << vector.y << ", " << vector.z
                << ")";
}

} // namespace gmt::hermite
