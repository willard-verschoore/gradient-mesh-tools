#include <QVector2D>
#include <QVector3D>
#include <cmath>

/// A tuple of (coordinates, color), acting as a 5D vector.
/**
 * A single element of a patch's control matrix,
 * consisting of a 2D position and an RGB color.
 */
struct Interpolant
{
  /// Initializes all components to zero.
  constexpr Interpolant() = default;
  /// Initializes all components to the given value.
  constexpr Interpolant(float x) : coords(x, x), color(x, x, x) {}
  /// Initializes coordinates to the given value, and color to zero.
  constexpr Interpolant(QVector2D coords) : coords(coords), color() {}
  /// Initializes coordinates and color to the given values.
  constexpr Interpolant(QVector2D coords, QVector3D color)
      : coords(coords), color(color)
  {
  }
  /// Initializes all five components to the specified values.
  constexpr Interpolant(float x, float y, float r, float g, float b)
      : coords(x, y), color(r, g, b)
  {
  }

  /// Number of components in this vector.
  static constexpr int COMPONENTS = 5;

  /// Index-based access of the vector's components.
  /**
   * Indices 0-1 represent coordinates,
   * whereas 2-4 represent color.
   */
  float& operator[](int index);
  /// Index-based access of the vector's components.
  /**
   * Indices 0-1 represent coordinates,
   * whereas 2-4 represent color.
   */
  float operator[](int index) const;

  /// The vector's coordinates in world space.
  QVector2D coords;
  /// The vector's color in RGB colorspace.
  QVector3D color;
};

/// Element-wise negates all components of a vector.
Interpolant operator-(Interpolant i);
/// Scales an interpolant by the inverse of the given value.
Interpolant& operator/=(Interpolant& left, float right);
/// Scales an interpolant by the given value.
Interpolant& operator*=(Interpolant& left, float right);
/// Element-wise multiplication between two interpolants.
Interpolant& operator*=(Interpolant& left, const Interpolant& right);
/// Element-wise addition between two interpolants.
Interpolant& operator+=(Interpolant& left, const Interpolant& right);
/// Element-wise subtraction between two interpolants.
Interpolant& operator-=(Interpolant& left, const Interpolant& right);

/// Scales an interpolant by the inverse of the given value.
inline Interpolant operator/(Interpolant left, float right)
{
  return left /= right;
}
/// Scales an interpolant by the given value.
inline Interpolant operator*(Interpolant left, float right)
{
  return left *= right;
}
/// Scales an interpolant by the given value.
inline Interpolant operator*(float left, Interpolant right)
{
  return right * left;
}
/// Element-wise multiplication between two interpolants.
inline Interpolant operator*(Interpolant left, const Interpolant& right)
{
  return left *= right;
}
/// Element-wise addition between two interpolants.
inline Interpolant operator+(Interpolant left, const Interpolant& right)
{
  return left += right;
}
/// Element-wise subtraction between two interpolants.
inline Interpolant operator-(Interpolant left, const Interpolant& right)
{
  return left -= right;
}
/// Fuzzy comparison of interpolant with fully zeroed interpolant.
bool is_zero(const Interpolant& i, float coords_epsilon, float color_epsilon);
/// The 'length' of the 5D vector of an Interpolant.
float length(const Interpolant& i);
