#pragma once

#include <ostream>

struct Vector2
{
  union
  {
    float data[2];

    struct
    {
      float x, y;
    };
  };

  constexpr Vector2(float x = 0.0f, float y = 0.0f) : x(x), y(y) {}

  Vector2 &add(Vector2 const &other);
  Vector2 &add(float value);

  Vector2 &operator+=(Vector2 const &other);
  Vector2 &operator+=(float value);

  Vector2 &subtract(Vector2 const &other);
  Vector2 &subtract(float value);

  Vector2 &operator-=(Vector2 const &other);
  Vector2 &operator-=(float value);

  Vector2 &multiply(Vector2 const &other);
  Vector2 &multiply(float value);

  Vector2 &operator*=(Vector2 const &value);
  Vector2 &operator*=(float value);

  Vector2 &divide(Vector2 const &other);
  Vector2 &divide(float value);

  Vector2 &operator/=(Vector2 const &value);
  Vector2 &operator/=(float value);

  float &operator[](int index);
  float operator[](int index) const;

  float length() const;
  float length_squared() const;

  Vector2 normalized() const;
  void normalize();
};

Vector2 add(Vector2 const &left, Vector2 const &right);
Vector2 add(Vector2 const &left, float right);
Vector2 add(float left, Vector2 const &right);

Vector2 operator+(Vector2 const &left, Vector2 const &right);
Vector2 operator+(Vector2 const &left, float right);
Vector2 operator+(float left, Vector2 const &right);

Vector2 subtract(Vector2 const &left, Vector2 const &right);
Vector2 subtract(Vector2 const &left, float right);
Vector2 subtract(float left, Vector2 const &right);

Vector2 operator-(Vector2 const &left, Vector2 const &right);
Vector2 operator-(Vector2 const &left, float right);
Vector2 operator-(float left, Vector2 const &right);

Vector2 multiply(Vector2 const &left, Vector2 const &right);
Vector2 multiply(Vector2 const &left, float right);
Vector2 multiply(float left, Vector2 const &right);

Vector2 operator*(Vector2 const &left, Vector2 const &right);
Vector2 operator*(Vector2 const &left, float right);
Vector2 operator*(float left, Vector2 const &right);

Vector2 divide(Vector2 const &left, Vector2 const &right);
Vector2 divide(Vector2 const &left, float right);
Vector2 divide(float left, Vector2 const &right);

Vector2 operator/(Vector2 const &left, Vector2 const &right);
Vector2 operator/(Vector2 const &left, float right);
Vector2 operator/(float left, Vector2 const &right);

Vector2 negate(Vector2 const &value);
Vector2 operator-(Vector2 const &value);

float dot(Vector2 const &left, Vector2 const &right);

float distance(Vector2 const &left, Vector2 const &right);
float distance_squared(Vector2 const &left, Vector2 const &right);

std::ostream &operator<<(std::ostream &stream, Vector2 const &vector);
