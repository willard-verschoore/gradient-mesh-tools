#pragma once

#include <ostream>

/// A vector in 3D space.
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

  /// Constructs a 3D vector from the specified components.
  /**
   * @param x,y,z The vector's components.
   */
  constexpr Vector3(float x = 0.0f, float y = 0.0f, float z = 0.0f)
      : x(x), y(y), z(z)
  {
  }

  /// Adds the given vector to this vector.
  /**
   * @param other The Vector3D to add to this vector.
   * @return A reference to this vector.
   */
  Vector3 &add(Vector3 const &other);

  /// Adds the given value to each component of this vector.
  /**
   * @param value The value to add to each component of this vector.
   * @return A reference to this vector.
   */
  Vector3 &add(float value);

  /// Adds the given vector to this vector.
  /**
   * @param other The Vector3D to add to this vector.
   * @return A reference to this vector.
   */
  Vector3 &operator+=(Vector3 const &other);

  /// Adds the given value to each component of this vector.
  /**
   * @param value The value to add to each component of this vector.
   * @return A reference to this vector.
   */
  Vector3 &operator+=(float value);

  /// Subtracts the given vector from this vector.
  /**
   * @param other The Vector3D to subtract from this vector.
   * @return A reference to this vector.
   */
  Vector3 &subtract(Vector3 const &other);

  /// Subtracts the given value from each component of this vector.
  /**
   * @param value The value to subtract from each component of this vector.
   * @return A reference to this vector.
   */
  Vector3 &subtract(float value);

  /// Subtracts the given vector from this vector.
  /**
   * @param other The Vector3D to subtract from this vector.
   * @return A reference to this vector.
   */
  Vector3 &operator-=(Vector3 const &other);

  /// Subtracts the given value from each component of this vector.
  /**
   * @param value The value to subtract from each component of this vector.
   * @return A reference to this vector.
   */
  Vector3 &operator-=(float value);

  /// Multiplies this vector with \c other component-wise.
  /**
   * Multiplies each component of this vector with the corresponding component
   * of \c other. Note that this is not the same as dot() or cross().
   *
   * @param other The Vector3D to multiply this one with.
   * @return A reference to this vector.
   */
  Vector3 &multiply(Vector3 const &other);

  /// Multiplies this vector with the given value.
  /**
   * @param value The value to multiply this vector with.
   * @return A reference to this vector.
   */
  Vector3 &multiply(float value);

  /// Multiplies this vector with \c other component-wise.
  /**
   * Multiplies each component of this vector with the corresponding component
   * of \c other. Note that this is not the same as dot() or cross().
   *
   * @param other The Vector3D to multiply this one with.
   * @return A reference to this vector.
   */
  Vector3 &operator*=(Vector3 const &value);

  /// Multiplies this vector with the given value.
  /**
   * @param value The value to multiply this vector with.
   * @return A reference to this vector.
   */
  Vector3 &operator*=(float value);

  /// Divides this vector by \c other component-wise.
  /**
   * Divides each component of this vector by the corresponding component of
   * \c other.
   *
   * @param other The Vector3D to divide this one by.
   * @return A reference to this vector.
   */
  Vector3 &divide(Vector3 const &other);

  /// Divides this vector by the given value.
  /**
   * @param value The value to divide this vector by.
   * @return A reference to this vector.
   */
  Vector3 &divide(float value);

  /// Divides this vector by \c other component-wise.
  /**
   * Divides each component of this vector by the corresponding component of
   * \c other.
   *
   * @param other The Vector3D to divide this one by.
   * @return A reference to this vector.
   */
  Vector3 &operator/=(Vector3 const &value);

  /// Divides this vector by the given value.
  /**
   * @param value The value to divide this vector by.
   * @return A reference to this vector.
   */
  Vector3 &operator/=(float value);

  /// Returns a reference to the component at the specified coordinate.
  /**
   * Note that \c index should be a value between 0 and 2 for a Vector3D. No
   * bounds checks are performed.
   *
   * @param index The index of the component to access.
   * @return A reference to the component at \c index.
   */
  float &operator[](int index);

  /// Returns a copy of the component at the specified coordinate.
  /**
   * Note that \c index should be a value between 0 and 2 for a Vector3D. No
   * bounds checks are performed.
   *
   * @param index The index of the component to access.
   * @return A copy of the component at \c index.
   */
  float operator[](int index) const;

  /// Returns the length of this vector.
  /**
   * @return The length of this vector.
   */
  float length() const;

  /// Returns the squared length of this vector.
  /**
   * This is faster to compute than length(), so it is preferred in cases where
   * only the relative length of vectors matters.
   *
   * Equivalent to applying dot() with this vector as both arguments.
   *
   * @return The squared length of this vector.
   */
  float length_squared() const;

  /// Returns a normalized copy of this vector.
  /**
   * @return A copy of this vector normalized to have length 1.
   */
  Vector3 normalized() const;

  /// Normalize this vector.
  void normalize();
};

/// Returns the addition of the \c left and \c right vectors.
/**
 * @param left,right The Vector3D objects to add together.
 * @return The result of adding \c left and \c right.
 */
Vector3 add(Vector3 const &left, Vector3 const &right);

/// Returns the addition of the \c left vector and the \c right value.
/**
 * @param left The Vector3D object to add.
 * @param right The value to add to each component of \c left.
 * @return The result of adding \c left and \c right.
 */
Vector3 add(Vector3 const &left, float right);

/// Returns the addition of the \c left value and the \c right vector.
/**
 * @param left The value to add to each component of \c right.
 * @param right The Vector3D object to add.
 * @return The result of adding \c left and \c right.
 */
Vector3 add(float left, Vector3 const &right);

/// Returns the addition of the \c left and \c right vectors.
/**
 * @param left,right The Vector3D objects to add together.
 * @return The result of adding \c left and \c right.
 */
Vector3 operator+(Vector3 const &left, Vector3 const &right);

/// Returns the addition of the \c left vector and the \c right value.
/**
 * @param left The Vector3D object to add.
 * @param right The value to add to each component of \c left.
 * @return The result of adding \c left and \c right.
 */
Vector3 operator+(Vector3 const &left, float right);

/// Returns the addition of the \c left value and the \c right vector.
/**
 * @param left The value to add to each component of \c right.
 * @param right The Vector3D object to add.
 * @return The result of adding \c left and \c right.
 */
Vector3 operator+(float left, Vector3 const &right);

/// Returns the result of subtracting \c right from \c left.
/**
 * @param left The Vector3D object on the left side of the subtraction.
 * @param right The Vector3D object on the right side of the subtraction.
 * @return The result of subtracting \c right from \c left.
 */
Vector3 subtract(Vector3 const &left, Vector3 const &right);

/// Returns the result of subtracting \c right from \c left.
/**
 * @param left The Vector3D object on the left side of the subtraction.
 * @param right The value to subtract from each component of \c left.
 * @return The result of subtracting \c right from \c left.
 */
Vector3 subtract(Vector3 const &left, float right);

// TODO: Remove, this is silly.
Vector3 subtract(float left, Vector3 const &right);

/// Returns the result of subtracting \c right from \c left.
/**
 * @param left The Vector3D object on the left side of the subtraction.
 * @param right The Vector3D object on the right side of the subtraction.
 * @return The result of subtracting \c right from \c left.
 */
Vector3 operator-(Vector3 const &left, Vector3 const &right);

/// Returns the result of subtracting \c right from \c left.
/**
 * @param left The Vector3D object on the left side of the subtraction.
 * @param right The value to subtract from each component of \c left.
 * @return The result of subtracting \c right from \c left.
 */
Vector3 operator-(Vector3 const &left, float right);

// TODO: Remove, this is silly.
Vector3 operator-(float left, Vector3 const &right);

/// Returns the component-wise multiplication of \c left and \c right.
/**
 * @param left,right The Vector3D objects to multiply together.
 * @return The component-wise multiplication of \c left and \c right.
 */
Vector3 multiply(Vector3 const &left, Vector3 const &right);

/// Returns the multiplication of the \c left vector with the \c right value.
/**
 * @param left The Vector3D object to be multiplied.
 * @param right The value to multiply each component of \c right with.
 * @return The result of multiplying \c left with \c right.
 */
Vector3 multiply(Vector3 const &left, float right);

/// Returns the multiplication of the \c right vector with the \c left value.
/**
 * @param left The value to multiply each component of \c left with.
 * @param right The Vector3D object to be multiplied.
 * @return The result of multiplying \c right with \c left.
 */
Vector3 multiply(float left, Vector3 const &right);

/// Returns the component-wise multiplication of \c left and \c right.
/**
 * @param left,right The Vector3D objects to multiply together.
 * @return The component-wise multiplication of \c left and \c right.
 */
Vector3 operator*(Vector3 const &left, Vector3 const &right);

/// Returns the multiplication of the \c left vector with the \c right value.
/**
 * @param left The Vector3D object to be multiplied.
 * @param right The value to multiply each component of \c right with.
 * @return The result of multiplying \c left with \c right.
 */
Vector3 operator*(Vector3 const &left, float right);

/// Returns the multiplication of the \c right vector with the \c left value.
/**
 * @param left The value to multiply each component of \c left with.
 * @param right The Vector3D object to be multiplied.
 * @return The result of multiplying \c right with \c left.
 */
Vector3 operator*(float left, Vector3 const &right);

/// Returns the component-wise division of \c left by \c right.
/**
 * @param left The Vector3D object on the left side of the subtraction.
 * @param right The Vector3D object on the right side of the subtraction.
 * @return The result of dividing \c left by \c right.
 */
Vector3 divide(Vector3 const &left, Vector3 const &right);

/// Returns the division of \c left by \c right.
/**
 * @param left The Vector3D object on the left side of the subtraction.
 * @param right The value to divide each component of \c left by.
 * @return The result of dividing \c left by \c right.
 */
Vector3 divide(Vector3 const &left, float right);

// TODO: Remove, this is silly.
Vector3 divide(float left, Vector3 const &right);

/// Returns the component-wise division of \c left by \c right.
/**
 * @param left The Vector3D object on the left side of the subtraction.
 * @param right The Vector3D object on the right side of the subtraction.
 * @return The result of dividing \c left by \c right.
 */
Vector3 operator/(Vector3 const &left, Vector3 const &right);

/// Returns the division of \c left by \c right.
/**
 * @param left The Vector3D object on the left side of the subtraction.
 * @param right The value to divide each component of \c left by.
 * @return The result of dividing \c left by \c right.
 */
Vector3 operator/(Vector3 const &left, float right);

// TODO: Remove, this is silly.
Vector3 operator/(float left, Vector3 const &right);

/// Returns the negation of the given vector.
/**
 * @param value The Vector3D object to negate.
 * @return The result of negating each component of \c value.
 */
Vector3 negate(Vector3 const &value);

/// Returns the negation of the given vector.
/**
 * @param value The Vector3D object to negate.
 * @return The result of negating each component of \c value.
 */
Vector3 operator-(Vector3 const &value);

/// Returns the dot product of \c left and \c right.
/**
 * @param left,right The Vector3D objects to take the dot product of.
 * @return The dot product of \c left and \c right.
 */
float dot(Vector3 const &left, Vector3 const &right);

/// Returns the dot product of \c left and \c right.
/**
 * @param left,right The Vector3D objects to take the cross product of.
 * @return The cross product of \c left and \c right.
 */
Vector3 cross(Vector3 const &left, Vector3 const &right);

/// Returns the distance between \c left and \c right.
/**
 * Equivalent to taking the length() of the difference.
 *
 * @param left,right The Vector3D objects get the distance between.
 * @return The distance between \c left and \c right.
 */
float distance(Vector3 const &left, Vector3 const &right);

/// Returns the squared distance between \c left and \c right.
/**
 * Equivalent to taking the length_squared() of the difference.
 *
 * @param left,right The Vector3D objects get the squared distance between.
 * @return The squared distance between \c left and \c right.
 */
float distance_squared(Vector3 const &left, Vector3 const &right);

/// Writes the given vector to the given stream.
/**
 * Will write the vector as (x, y, z), i.e., a space separated list of the
 * vector's components between parentheses.
 *
 * @param stream The output stream to write to.
 * @param vector The Vector3D object to write.
 * @return A reference to \c stream.
 */
std::ostream &operator<<(std::ostream &stream, Vector3 const &vector);
