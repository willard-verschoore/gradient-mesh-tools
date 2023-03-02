#pragma once

#include <algorithm>

namespace gmt
{

/// An ordered pair of (start, end) floating-point values.
/**
 * An ordered pair of (start, end) floating-point values,
 * representing an interval of real numbers.
 */
struct Interval
{
  Interval(float start, float end) : start(start), end(end) {}

  /// Linearly interpolates the two endpoints of the interval at point t.
  /**
   * @param t: the point along the interval. Must be between 0 and 1.
   */
  float operator()(float t) const { return start + t * (end - start); };

  float start, end;
};

/// The extent (or length) of the given interval.
inline float length(const Interval& i) { return i.end - i.start; }

/// Whether or not the interval has an extent equal to 0.
inline bool is_empty(const Interval& i) { return i.end == i.start; }

inline bool is_full(const Interval& i) { return length(i) > 0.99; }

/// Whether point `t` is contained in the given interval.
inline bool contains(const Interval& i, float t)
{
  return t >= i.start && t < i.end;
}

/// Returns true if the two intervals overlap, or false if they're disjoint.
inline bool intersects(const Interval& a, const Interval& b)
{
  return a.start <= b.end && a.end >= b.start;
}

/**
 * Returns the equivalent interval when flipping the number line,
 * thus mapping [-inf, inf] to [inf, -inf].
 */
inline Interval converse(Interval i)
{
  std::swap(i.start, i.end);
  i.start = 1.0f - i.start;
  i.end = 1.0f - i.end;
  return i;
}

/**
 * Returns the equivalent interval when mapping its containing interval
 * from [0, 1] to `new_outer`.
 */
inline Interval reparametrize(Interval inner, const Interval& new_outer)
{
  inner.start = (inner.start - new_outer.start) / length(new_outer);
  inner.end = (inner.end - new_outer.start) / length(new_outer);
  return inner;
}

/**
 * Returns the equivalent position t when mapping the parent_t to the child
 * interval.
 */
inline float relative_child_position(Interval child, float parent_t)
{
  return (parent_t - child.start) / length(child);
}

} // namespace gmt
