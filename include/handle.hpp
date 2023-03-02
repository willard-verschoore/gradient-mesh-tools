#pragma once

#include <variant>

#include "hermite/hermite.hpp"
#include "storage.hpp"

struct HalfEdge;

/// A user-editable tangent handle, used to edit the shape of a curve.
/**
 * In the UI, tangent handles are shown in Bezier form,
 * but they're internally represented in hermite form.
 * Thus, when dragging, the drag amount is first scaled to match.
 */
struct Handle
{
  Handle(hermite::Interpolant tangent) : tangent(tangent) {}

  /// Tangent of the curve in Hermite form.
  /**
   * This tangent does not necessarily correspond to the patch's
   * tangent vectors in its control matrix, but may be flipped
   * in the opposite direction.
   * This is to make sure that opposite half-edges are
   * completely symmetric between each other.
   */
  hermite::Interpolant tangent;
  /// Back-reference to the edge owning this handle.
  Id<HalfEdge> edge;
};
