#pragma once

#include <ostream>

namespace gmt::hermite
{

/// A vector in 2D space.
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

  /// Constructs a 2D vector from the specified components.
  /**
   * @param x,y The vector's components.
   */
  constexpr Vector2(float x = 0.0f, float y = 0.0f) : x(x), y(y) {}

  /// Adds the given vector to this vector.
  /**
   * @param other The Vector2 to add to this vector.
   * @return A reference to this vector.
   */
  Vector2 &add(Vector2 const &other);

  /// Adds the given value to each component of this vector.
  /**
   * @param value The value to add to each component of this vector.
   * @return A reference to this vector.
   */
  Vector2 &add(float value);

  /// Adds the given vector to this vector.
  /**
   * @param other The Vector2 to add to this vector.
   * @return A reference to this vector.
   */
  Vector2 &operator+=(Vector2 const &other);

  /// Adds the given value to each component of this vector.
  /**
   * @param value The value to add to each component of this vector.
   * @return A reference to this vector.
   */
  Vector2 &operator+=(float value);

  /// Subtracts the given vector from this vector.
  /**
   * @param other The Vector2 to subtract from this vector.
   * @return A reference to this vector.
   */
  Vector2 &subtract(Vector2 const &other);

  /// Subtracts the given value from each component of this vector.
  /**
   * @param value The value to subtract from each component of this vector.
   * @return A reference to this vector.
   */
  Vector2 &subtract(float value);

  /// Subtracts the given vector from this vector.
  /**
   * @param other The Vector2 to subtract from this vector.
   * @return A reference to this vector.
   */
  Vector2 &operator-=(Vector2 const &other);

  /// Subtracts the given value from each component of this vector.
  /**
   * @param value The value to subtract from each component of this vector.
   * @return A reference to this vector.
   */
  Vector2 &operator-=(float value);

  /// Multiplies this vector with \c other component-wise.
  /**
   * Multiplies each component of this vector with the corresponding component
   * of \c other. Note that this is not the same as dot().
   *
   * @param other The Vector2 to multiply this one with.
   * @return A reference to this vector.
   */
  Vector2 &multiply(Vector2 const &other);

  /// Multiplies this vector with the given value.
  /**
   * @param value The value to multiply this vector with.
   * @return A reference to this vector.
   */
  Vector2 &multiply(float value);

  /// Multiplies this vector with \c other component-wise.
  /**
   * Multiplies each component of this vector with the corresponding component
   * of \c other. Note that this is not the same as dot().
   *
   * @param other The Vector2 to multiply this one with.
   * @return A reference to this vector.
   */
  Vector2 &operator*=(Vector2 const &value);

  /// Multiplies this vector with the given value.
  /**
   * @param value The value to multiply this vector with.
   * @return A reference to this vector.
   */
  Vector2 &operator*=(float value);

  /// Divides this vector by \c other component-wise.
  /**
   * Divides each component of this vector by the corresponding component of
   * \c other.
   *
   * @param other The Vector2 to divide this one by.
   * @return A reference to this vector.
   */
  Vector2 &divide(Vector2 const &other);

  /// Divides this vector by the given value.
  /**
   * @param value The value to divide this vector by.
   * @return A reference to this vector.
   */
  Vector2 &divide(float value);

  /// Divides this vector by \c other component-wise.
  /**
   * Divides each component of this vector by the corresponding component of
   * \c other.
   *
   * @param other The Vector2 to divide this one by.
   * @return A reference to this vector.
   */
  Vector2 &operator/=(Vector2 const &value);

  /// Divides this vector by the given value.
  /**
   * @param value The value to divide this vector by.
   * @return A reference to this vector.
   */
  Vector2 &operator/=(float value);

  /// Returns a reference to the component at the specified coordinate.
  /**
   * Note that \c index should be a value between 0 and 1 for a Vector2. No
   * bounds checks are performed.
   *
   * @param index The index of the component to access.
   * @return A reference to the component at \c index.
   */
  float &operator[](int index);

  /// Returns a copy of the component at the specified coordinate.
  /**
   * Note that \c index should be a value between 0 and 1 for a Vector2. No
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
  Vector2 normalized() const;

  /// Normalize this vector.
  void normalize();
};

/// Returns the addition of the \c left and \c right vectors.
/**
 * @param left,right The Vector2 objects to add together.
 * @return The result of adding \c left and \c right.
 */
Vector2 add(Vector2 const &left, Vector2 const &right);

/// Returns the addition of the \c left vector and the \c right value.
/**
 * @param left The Vector2 object to add.
 * @param right The value to add to each component of \c left.
 * @return The result of adding \c left and \c right.
 */
Vector2 add(Vector2 const &left, float right);

/// Returns the addition of the \c left value and the \c right vector.
/**
 * @param left The value to add to each component of \c right.
 * @param right The Vector2 object to add.
 * @return The result of adding \c left and \c right.
 */
Vector2 add(float left, Vector2 const &right);

/// Returns the addition of the \c left and \c right vectors.
/**
 * @param left,right The Vector2 objects to add together.
 * @return The result of adding \c left and \c right.
 */
Vector2 operator+(Vector2 const &left, Vector2 const &right);

/// Returns the addition of the \c left vector and the \c right value.
/**
 * @param left The Vector2 object to add.
 * @param right The value to add to each component of \c left.
 * @return The result of adding \c left and \c right.
 */
Vector2 operator+(Vector2 const &left, float right);

/// Returns the addition of the \c left value and the \c right vector.
/**
 * @param left The value to add to each component of \c right.
 * @param right The Vector2 object to add.
 * @return The result of adding \c left and \c right.
 */
Vector2 operator+(float left, Vector2 const &right);

/// Returns the result of subtracting \c right from \c left.
/**
 * @param left The Vector2 object on the left side of the subtraction.
 * @param right The Vector2 object on the right side of the subtraction.
 * @return The result of subtracting \c right from \c left.
 */
Vector2 subtract(Vector2 const &left, Vector2 const &right);

/// Returns the result of subtracting \c right from \c left.
/**
 * @param left The Vector2 object on the left side of the subtraction.
 * @param right The value to subtract from each component of \c left.
 * @return The result of subtracting \c right from \c left.
 */
Vector2 subtract(Vector2 const &left, float right);

/// Returns the result of subtracting \c right from \c left.
/**
 * @param left The Vector2 object on the left side of the subtraction.
 * @param right The Vector2 object on the right side of the subtraction.
 * @return The result of subtracting \c right from \c left.
 */
Vector2 operator-(Vector2 const &left, Vector2 const &right);

/// Returns the result of subtracting \c right from \c left.
/**
 * @param left The Vector2 object on the left side of the subtraction.
 * @param right The value to subtract from each component of \c left.
 * @return The result of subtracting \c right from \c left.
 */
Vector2 operator-(Vector2 const &left, float right);

/// Returns the component-wise multiplication of \c left and \c right.
/**
 * @param left,right The Vector2 objects to multiply together.
 * @return The component-wise multiplication of \c left and \c right.
 */
Vector2 multiply(Vector2 const &left, Vector2 const &right);

/// Returns the multiplication of the \c left vector with the \c right value.
/**
 * @param left The Vector2 object to be multiplied.
 * @param right The value to multiply each component of \c right with.
 * @return The result of multiplying \c left with \c right.
 */
Vector2 multiply(Vector2 const &left, float right);

/// Returns the multiplication of the \c right vector with the \c left value.
/**
 * @param left The value to multiply each component of \c left with.
 * @param right The Vector2 object to be multiplied.
 * @return The result of multiplying \c right with \c left.
 */
Vector2 multiply(float left, Vector2 const &right);

/// Returns the component-wise multiplication of \c left and \c right.
/**
 * @param left,right The Vector2 objects to multiply together.
 * @return The component-wise multiplication of \c left and \c right.
 */
Vector2 operator*(Vector2 const &left, Vector2 const &right);

/// Returns the multiplication of the \c left vector with the \c right value.
/**
 * @param left The Vector2 object to be multiplied.
 * @param right The value to multiply each component of \c right with.
 * @return The result of multiplying \c left with \c right.
 */
Vector2 operator*(Vector2 const &left, float right);

/// Returns the multiplication of the \c right vector with the \c left value.
/**
 * @param left The value to multiply each component of \c left with.
 * @param right The Vector2 object to be multiplied.
 * @return The result of multiplying \c right with \c left.
 */
Vector2 operator*(float left, Vector2 const &right);

/// Returns the component-wise division of \c left by \c right.
/**
 * @param left The Vector2 object on the left side of the subtraction.
 * @param right The Vector2 object on the right side of the subtraction.
 * @return The result of dividing \c left by \c right.
 */
Vector2 divide(Vector2 const &left, Vector2 const &right);

/// Returns the division of \c left by \c right.
/**
 * @param left The Vector2 object on the left side of the subtraction.
 * @param right The value to divide each component of \c left by.
 * @return The result of dividing \c left by \c right.
 */
Vector2 divide(Vector2 const &left, float right);

/// Returns the component-wise division of \c left by \c right.
/**
 * @param left The Vector2 object on the left side of the subtraction.
 * @param right The Vector2 object on the right side of the subtraction.
 * @return The result of dividing \c left by \c right.
 */
Vector2 operator/(Vector2 const &left, Vector2 const &right);

/// Returns the division of \c left by \c right.
/**
 * @param left The Vector2 object on the left side of the subtraction.
 * @param right The value to divide each component of \c left by.
 * @return The result of dividing \c left by \c right.
 */
Vector2 operator/(Vector2 const &left, float right);

/// Returns the negation of the given vector.
/**
 * @param value The Vector2 object to negate.
 * @return The result of negating each component of \c value.
 */
Vector2 negate(Vector2 const &value);

/// Returns the negation of the given vector.
/**
 * @param value The Vector2 object to negate.
 * @return The result of negating each component of \c value.
 */
Vector2 operator-(Vector2 const &value);

/// Returns the dot product of \c left and \c right.
/**
 * @param left,right The Vector2 objects to take the dot product of.
 * @return The dot product of \c left and \c right.
 */
float dot(Vector2 const &left, Vector2 const &right);

/// Returns the distance between \c left and \c right.
/**
 * Equivalent to taking the length() of the difference.
 *
 * @param left,right The Vector2 objects to get the distance between.
 * @return The distance between \c left and \c right.
 */
float distance(Vector2 const &left, Vector2 const &right);

/// Returns the squared distance between \c left and \c right.
/**
 * Equivalent to taking the length_squared() of the difference.
 *
 * @param left,right The Vector2 objects to get the squared distance between.
 * @return The squared distance between \c left and \c right.
 */
float distance_squared(Vector2 const &left, Vector2 const &right);

/// Checks whether \c left and \c right are equal.
/**
 * Two Vector2 objects are equal if each corresponding pair of components is
 * equal.
 *
 * @param left,right The Vector2 objects to compare.
 * @return Whether \c left and \c right are equal.
 */
bool operator==(Vector2 const &left, Vector2 const &right);

/// Checks whether \c left and \c right are not equal.
/**
 * Two Vector2 objects are equal if each corresponding pair of components is
 * equal. They are therefore unequal if there is at least one unequal pair.
 *
 * @param left,right The Vector2 objects to compare.
 * @return Whether \c left and \c right are not equal.
 */
bool operator!=(Vector2 const &left, Vector2 const &right);

/// Writes the given vector to the given stream.
/**
 * Will write the vector as (x, y), i.e., a space separated list of the vector's
 * components between parentheses.
 *
 * @param stream The output stream to write to.
 * @param vector The Vector2 object to write.
 * @return A reference to \c stream.
 */
std::ostream &operator<<(std::ostream &stream, Vector2 const &vector);

} // namespace gmt::hermite
