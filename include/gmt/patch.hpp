#pragma once

#include "storage.hpp"

namespace gmt
{

struct HalfEdge;

/// One individual bicubic hermite patch in the gradient mesh.
struct Patch
{
  /// Reference to one of the edges bordering the patch.
  Id<HalfEdge> side;
};

} // namespace gmt
