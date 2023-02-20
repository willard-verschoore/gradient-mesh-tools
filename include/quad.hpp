#pragma once

#include <functional>
#include <utility>

/// A generic class storing exactly four instances of T.
/**
 * A generic container class that stores four values of type T,
 * one for each side of a quad.
 */
template <typename T>
struct Quad
{
  Quad(T top, T left, T bottom, T right)
      : top(std::move(top)),
        left(std::move(left)),
        bottom(std::move(bottom)),
        right(std::move(right))
  {
  }

  T top, left, bottom, right;
};

/// Applies a function to all elements of the quad, returning the result.
template <typename T, typename U>
Quad<U> map(const std::function<U(const T&)>& fn, const Quad<T>& q)
{
  return Quad<U>(fn(q.top), fn(q.left), fn(q.bottom), fn(q.right));
}
