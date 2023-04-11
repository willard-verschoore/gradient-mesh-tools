#pragma once

#include <map>
#include <optional>
#include <variant>

#include "control-point.hpp"
#include "handle.hpp"
#include "hermite/hermite.hpp"
#include "interval.hpp"
#include "patch.hpp"

namespace gmt
{

// Forward declarations
struct HalfEdge;

/// A half-edge that is connected to a control point.
/**
 * A "Parent" half-edge, which is directly
 * connected to a user-editable control point.
 */
struct Parent
{
  /// Control point (vertex) that this half edge originates from.
  Id<ControlPoint> origin;
  /// Tangent of the curve at the start and end points.
  std::array<Id<Handle>, 2> handles;
};

/// A half-edge that computes its parameters from a parent.
/**
 * A child half-edge, whose geometry is computed
 * from its parent half-edge.
 */
struct Child
{
  /// The parent Half-edge from which data is computed.
  Id<HalfEdge> parent;

  /// The interval along the curve on which the edge's origin lies.
  /**
   * If this interval is empty (start and endpoint are the same)
   * the curve is orthogonal to its parent.
   */
  Interval interval;

  /// Tangent of the curve at the start and end points.
  /**
   * If not present (std::nullopt) this data is instead computed
   */
  std::optional<std::array<Id<Handle>, 2>> handles;

  /// Whether this child has been given tangent handles for recoloring.
  /**
   * In order to properly recolor child edges the color of its tangents needs to
   * be set as well. For a child without tangent handles this means that new
   * handles need to be generated. This flag signals that the child's handles
   * have been generated for this purpose and should therefore only affect the
   * color of the edge while its position is still determined by its parent.
   */
  bool recolored = false;
};

/// An individual, augmented half-edge for use with `GradientMesh`es.
struct HalfEdge
{
  /// Creates a parent half-edge with the given parameters.
  HalfEdge(Id<ControlPoint> origin, std::array<Id<Handle>, 2> handles,
           hermite::Vector3 color,
           hermite::Interpolant twist = hermite::Interpolant(0.0f),
           std::optional<Id<HalfEdge>> twin = std::nullopt)
      : kind(Parent{origin, handles}), color(color), twist(twist), twin(twin)
  {
  }

  /// Creates a child half-edge with the given parameters.
  HalfEdge(Id<HalfEdge> parent, Interval interval,
           std::optional<std::array<Id<Handle>, 2>> handles,
           hermite::Vector3 color,
           hermite::Interpolant twist = hermite::Interpolant(0.0f),
           std::optional<Id<HalfEdge>> twin = std::nullopt)
      : kind(Child{parent, interval, handles}),
        color(color),
        twist(twist),
        twin(twin)
  {
  }

  /// The edge's interval relative to its parent, or [0, 1] if it has none.
  Interval interval() const;
  /// Same as interval(), except it also returns [0, 1] if the half-edge is
  /// orthogonal to its parent.
  Interval relative_interval() const;

  /// Retrieves the half-edge's tangent handles, in the order `[start, end]`
  std::optional<std::array<Id<Handle>, 2>> handles() const;
  /// Sets the half-edge's tangent handles.
  /**
   * @param handles: the two new handles, in the order `[start, end]`
   */
  void set_handles(const std::array<Id<Handle>, 2>& handles);

  /**
   * Returns whether splitting this half-edge will result
   * in a new level of the half-edge tree.
   */
  bool is_parentable() const;

  int level;
  /// Back-reference to the face this half edge is part of.
  /**
   * Typically set after construction.
   *  Always remember to initialize it with the `connect`
   *  function!
   */
  Id<Patch> patch;
  /// Next and previous half-edges in the face.
  /**
   * Typically set after construction.
   *  Always remember to initialize it with the `connect`
   *  function!
   */
  Id<HalfEdge> next, prev;
  /// Distinguishes whether this half-edge is a parent or a child.
  std::variant<Parent, Child> kind;
  /**
   * The color at the origin, stored in the half-edge
   * to enable sharp transitions.
   */
  hermite::Vector3 color;
  /// Twist vector at the origin point.
  /**
   * Should normally be 0,
   * but can be changed by splitting.
   */
  hermite::Interpolant twist = hermite::Interpolant(0.0f);
  /// Twin half-edge in the adjacent patch.
  /**
   * A half-edge may or may not have an opposite (twin) edge,
   * hence why the std::optional.
   */
  std::optional<Id<HalfEdge>> twin = std::nullopt;
  /**
   * If this edge has any children, this reference
   * points to its leftmost one (as in, the one whose
   * origin is the closest to the parent's).
   */
  std::optional<Id<HalfEdge>> leftmost_child = std::nullopt;
};

/// Sets a's "twin" edge to b, and b's "twin" edge to a.
void connect_twins(Id<HalfEdge> a, Id<HalfEdge> b, Storage<HalfEdge>& edges);
/// Sets edge's patch to 'patch', and patch's side to 'edge'.
void connect(Id<HalfEdge> edge, Id<Patch> patch, Storage<HalfEdge>& edges,
             Storage<Patch>& patches);

/**
 * Returns the edge opposite to the given one in its patch,
 * also known as its second-order neighbor.
 */
Id<HalfEdge> opposite(Id<HalfEdge> edge, const Storage<HalfEdge>& edges);

/// Returns the edge adjacent to the given one in the edge's outgoing direction.
/*
 * Computed as next(twin(next(edge))).
 */
std::optional<Id<HalfEdge>> adjacent_next(Id<HalfEdge> edge,
                                          const Storage<HalfEdge>& edges);

/// Returns the edge adjacent to the given one in the edge's reverse direction.
/*
 * Computed as prev(twin(prev(edge))).
 */
std::optional<Id<HalfEdge>> adjacent_prev(Id<HalfEdge> edge,
                                          const Storage<HalfEdge>& edges);

/** Returns the edge's nearest ancestor in the tree that a newly-split edge
 * should be parented to. In layman's terms, it returns the edge itself if it's
 * either a parent or orthogonal to its parent, or the edge's parent otherwise.
 * This is so that adjacent half-edges can occupy the same level in the edge
 * tree.
 */
Id<HalfEdge> eligible_parent(Id<HalfEdge> edge, const Storage<HalfEdge>& edges);

using ChildMap = std::map<float, Id<HalfEdge>>;
using ChildIterator = typename ChildMap::const_iterator;

/**
 * Finds all the children of the given edge,
 * returning them in an (ordered) map by their position relative to edge.
 * Only the parallel children are returned, not the orthogonal children.
 */
ChildMap children(Id<HalfEdge> edge, const Storage<HalfEdge>& edges);
/**
 * Finds all the siblings (children of its parent) of the given edge,
 * returning them in an (ordered) map by their position relative to edge.
 */
ChildMap siblings(Id<HalfEdge> edge, const Storage<HalfEdge>& edges);

/**
 * Finds all of the edges abutting the given one, which should be one more
 * than the number of T-junctions alongside the twin half-edge.
 * When is_relative is set to true, the relative position to the child
 * interval is used.
 */
ChildMap abutting_edges(Id<HalfEdge> edge, const Storage<HalfEdge>& edges,
                        bool is_relative = false);

/**
 * Returns the number of "sibling" half-edges
 * located alongside the neighboring twin half-edge.
 */
std::size_t siblings_abutting_twin(Id<HalfEdge> edge,
                                   const Storage<HalfEdge>& edges);

/**
 * Returns an iterator in children starting from child.
 * Starts from children.end() if child does not exist in
 * children.
 */
ChildIterator find(Id<HalfEdge> child, ChildMap& children);

/**
 * Returns the first (leftmost) child of an edge, or itself if it does not have
 * children.
 */
Id<HalfEdge> first_child(Id<HalfEdge> edge, const Storage<HalfEdge>& edges);

/**
 * Returns the last (rightmost) child of an edge, or itself if it does not have
 * children.
 */
Id<HalfEdge> last_child(Id<HalfEdge> edge, const Storage<HalfEdge>& edges);

size_t num_children(Id<HalfEdge> edge, const Storage<HalfEdge>& edges);

/// \cond HIDDEN_SYMBOLS
// Helper class for visiting a std::variant.
// Taken from https://en.cppreference.com/w/cpp/utility/variant/visit
template <class... Ts>
struct overloaded : Ts...
{
  using Ts::operator()...;
};
// Explicit deduction guide (to avoid specifying types)
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
/// \endcond

/**
 * Applies one of the two given functions to the edge's data, depending
 * on its kind, and returns its result. The two functions need to have
 * the same return type.
 * @param parent_fn: function to apply if the half-edge is a parent.
 * @param child_fn: function to apply if the half-edge is a child.
 * @param e: the half-edge to apply the functions to.
 */
template <typename F1, typename F2>
auto visit(F1&& parent_fn, F2&& child_fn, HalfEdge& e)
{
  return std::visit(
      overloaded{
          parent_fn,
          child_fn,
      },
      e.kind);
}

/**
 * Applies one of the two given functions to the edge's data, depending
 * on its kind, and returns its result. The two functions need to have
 * the same return type.
 * @param parent_fn: function to apply if the half-edge is a parent.
 * @param child_fn: function to apply if the half-edge is a child.
 * @param e: the half-edge to apply the functions to.
 */
template <typename F1, typename F2>
auto visit(F1&& parent_fn, F2&& child_fn, const HalfEdge& e)
{
  return std::visit(
      overloaded{
          parent_fn,
          child_fn,
      },
      e.kind);
}

} // namespace gmt
