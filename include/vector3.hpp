#pragma once

#include <ostream>

struct Vector3
{
  union
  {
    float data[3];

    struct
    {
      float x, y, z;
    };

    struct
    {
      float r, g, b;
    };
  };

  constexpr Vector3(float x = 0.0f, float y = 0.0f, float z = 0.0f)
      : x(x), y(y), z(z)
  {
  }

  Vector3 &add(Vector3 const &other);
  Vector3 &add(float value);

  Vector3 &operator+=(Vector3 const &other);
  Vector3 &operator+=(float value);

  Vector3 &subtract(Vector3 const &other);
  Vector3 &subtract(float value);

  Vector3 &operator-=(Vector3 const &other);
  Vector3 &operator-=(float value);

  Vector3 &multiply(Vector3 const &other);
  Vector3 &multiply(float value);

  Vector3 &operator*=(Vector3 const &value);
  Vector3 &operator*=(float value);

  Vector3 &divide(Vector3 const &other);
  Vector3 &divide(float value);

  Vector3 &operator/=(Vector3 const &value);
  Vector3 &operator/=(float value);

  float &operator[](int index);
  float operator[](int index) const;

  float length() const;
  float length_squared() const;

  Vector3 normalized() const;
  void normalize();
};

Vector3 add(Vector3 const &left, Vector3 const &right);
Vector3 add(Vector3 const &left, float right);
Vector3 add(float left, Vector3 const &right);

Vector3 operator+(Vector3 const &left, Vector3 const &right);
Vector3 operator+(Vector3 const &left, float right);
Vector3 operator+(float left, Vector3 const &right);

Vector3 subtract(Vector3 const &left, Vector3 const &right);
Vector3 subtract(Vector3 const &left, float right);
Vector3 subtract(float left, Vector3 const &right);

Vector3 operator-(Vector3 const &left, Vector3 const &right);
Vector3 operator-(Vector3 const &left, float right);
Vector3 operator-(float left, Vector3 const &right);

Vector3 multiply(Vector3 const &left, Vector3 const &right);
Vector3 multiply(Vector3 const &left, float right);
Vector3 multiply(float left, Vector3 const &right);

Vector3 operator*(Vector3 const &left, Vector3 const &right);
Vector3 operator*(Vector3 const &left, float right);
Vector3 operator*(float left, Vector3 const &right);

Vector3 divide(Vector3 const &left, Vector3 const &right);
Vector3 divide(Vector3 const &left, float right);
Vector3 divide(float left, Vector3 const &right);

Vector3 operator/(Vector3 const &left, Vector3 const &right);
Vector3 operator/(Vector3 const &left, float right);
Vector3 operator/(float left, Vector3 const &right);

Vector3 negate(Vector3 const &value);
Vector3 operator-(Vector3 const &value);

float dot(Vector3 const &left, Vector3 const &right);
Vector3 cross(Vector3 const &left, Vector3 const &right);

float distance(Vector3 const &left, Vector3 const &right);
float distance_squared(Vector3 const &left, Vector3 const &right);

std::ostream &operator<<(std::ostream &stream, Vector3 const &vector);
