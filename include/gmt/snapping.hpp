#pragma once

#include <utility>

namespace gmt
{

/// Default threshold when finding a parameter value in a std::map<float, T>.
inline float DEFAULT_SNAPPINESS = 0.01;

/// Fuzzy-finds an element in a std::map<float, T> or std::set<float>
/**
 * Works on both `const` and non-`const` associative containers.
 * @param container: the container to search.
 * @param target: the parameter value to search for.
 * @param epsilon: threshold to use when finding a value.
 * Taken from https://stackoverflow.com/a/13369134
 */
template <class Container>
auto find_approx(Container&& container, float target,
                 float epsilon = DEFAULT_SNAPPINESS)
    -> decltype(container.find(target))
{
  auto lower = container.lower_bound(target - epsilon);
  auto upper = container.upper_bound(target + epsilon);
  if (lower == upper)
  {
    // Item not found
    return container.end();
  }
  else
  {
    return lower;
  }
}

} // namespace gmt
