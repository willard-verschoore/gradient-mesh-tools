#include "vector2.hpp"

#include <cmath>

Vector2 &Vector2::add(Vector2 const &other)
{
  x += other.x;
  y += other.y;

  return *this;
}

Vector2 &Vector2::add(float value)
{
  x += value;
  y += value;

  return *this;
}

Vector2 &Vector2::operator+=(Vector2 const &other) { return add(other); }

Vector2 &Vector2::operator+=(float value) { return add(value); }

Vector2 add(Vector2 const &left, Vector2 const &right)
{
  Vector2 copy{left};
  return copy.add(right);
}

Vector2 add(Vector2 const &left, float right)
{
  Vector2 copy{left};
  return copy.add(right);
}

Vector2 add(float left, Vector2 const &right)
{
  Vector2 copy{right};
  return copy.add(left);
}

Vector2 operator+(Vector2 const &left, Vector2 const &right)
{
  return add(left, right);
}

Vector2 operator+(Vector2 const &left, float right) { return add(left, right); }

Vector2 operator+(float left, Vector2 const &right) { return add(left, right); }

Vector2 &Vector2::subtract(Vector2 const &other)
{
  x -= other.x;
  y -= other.y;

  return *this;
}

Vector2 &Vector2::subtract(float value)
{
  x -= value;
  y -= value;

  return *this;
}

Vector2 &Vector2::operator-=(Vector2 const &other) { return subtract(other); }

Vector2 &Vector2::operator-=(float value) { return subtract(value); }

Vector2 subtract(Vector2 const &left, Vector2 const &right)
{
  Vector2 copy{left};
  return copy.subtract(right);
}

Vector2 subtract(Vector2 const &left, float right)
{
  Vector2 copy{left};
  return copy.subtract(right);
}

Vector2 subtract(float left, Vector2 const &right)
{
  Vector2 copy{right};
  return copy.subtract(left);
}

Vector2 operator-(Vector2 const &left, Vector2 const &right)
{
  return subtract(left, right);
}

Vector2 operator-(Vector2 const &left, float right)
{
  return subtract(left, right);
}

Vector2 operator-(float left, Vector2 const &right)
{
  return subtract(left, right);
}

Vector2 &Vector2::multiply(Vector2 const &other)
{
  x *= other.x;
  y *= other.y;

  return *this;
}

Vector2 &Vector2::multiply(float value)
{
  x *= value;
  y *= value;

  return *this;
}

Vector2 &Vector2::operator*=(Vector2 const &other) { return multiply(other); }

Vector2 &Vector2::operator*=(float value) { return multiply(value); }

Vector2 multiply(Vector2 const &left, Vector2 const &right)
{
  Vector2 copy{left};
  return copy.multiply(right);
}

Vector2 multiply(Vector2 const &left, float right)
{
  Vector2 copy{left};
  return copy.multiply(right);
}

Vector2 multiply(float left, Vector2 const &right)
{
  Vector2 copy{right};
  return copy.multiply(left);
}

Vector2 operator*(Vector2 const &left, Vector2 const &right)
{
  return multiply(left, right);
}

Vector2 operator*(Vector2 const &left, float right)
{
  return multiply(left, right);
}

Vector2 operator*(float left, Vector2 const &right)
{
  return multiply(left, right);
}

Vector2 &Vector2::divide(Vector2 const &other)
{
  x /= other.x;
  y /= other.y;

  return *this;
}

Vector2 &Vector2::divide(float value)
{
  float inverse = 1.0f / value;

  x *= inverse;
  y *= inverse;

  return *this;
}

Vector2 &Vector2::operator/=(Vector2 const &other) { return divide(other); }

Vector2 &Vector2::operator/=(float value) { return divide(value); }

Vector2 divide(Vector2 const &left, Vector2 const &right)
{
  Vector2 copy{left};
  return copy.divide(right);
}

Vector2 divide(Vector2 const &left, float right)
{
  Vector2 copy{left};
  return copy.divide(right);
}

Vector2 divide(float left, Vector2 const &right)
{
  Vector2 copy{right};
  return copy.divide(left);
}

Vector2 operator/(Vector2 const &left, Vector2 const &right)
{
  return divide(left, right);
}

Vector2 operator/(Vector2 const &left, float right)
{
  return divide(left, right);
}

Vector2 operator/(float left, Vector2 const &right)
{
  return divide(left, right);
}

float &Vector2::operator[](int index) { return data[index]; }

float Vector2::operator[](int index) const { return data[index]; }

Vector2 negate(Vector2 const &value)
{
  Vector2 copy{value};

  copy.x = -copy.x;
  copy.y = -copy.y;

  return copy;
}

Vector2 operator-(Vector2 const &value) { return negate(value); }

float Vector2::length() const { return sqrtf(length_squared()); }

float Vector2::length_squared() const { return x * x + y * y; }

Vector2 Vector2::normalized() const { return (*this) / length(); }

void Vector2::normalize()
{
  float inverse = 1.0f / length();

  x *= inverse;
  y *= inverse;
}

float dot(Vector2 const &left, Vector2 const &right)
{
  return left.x * right.x + left.y * right.y;
}

float distance(Vector2 const &left, Vector2 const &right)
{
  return (left - right).length();
}

float distance_squared(Vector2 const &left, Vector2 const &right)
{
  return (left - right).length_squared();
}

std::ostream &operator<<(std::ostream &stream, Vector2 const &vector)
{
  return stream << "(" << vector.x << ", " << vector.y << ")";
}
