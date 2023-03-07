#pragma once

#include <array>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "half-edge.hpp"
#include "merge-options.hpp"
#include "tessellation.hpp"

namespace gmt
{

inline float MAXIMUM_TWIN_SNAPPING_OFFSET = 0.001;

/// The main Gradient Mesh data structure, represented as a half-edge graph.
/**
 * Functions related to splitting are implemented in the file `splitting.cpp`.
 * All others are implemented in `gradmesh.cpp`.
 */
class GradientMesh
{
 public:
  GradientMesh() = default;
  /// Create a basic gradient mesh.
  /**
   * The mesh will consist of a square centered at the origin.
   * @param side: The length of the squre's side.
   * @param colors: The color of the square's four corners,
   *    in the order (top left, top right, bottom right, bottom left).
   */
  GradientMesh(float side, const std::array<hermite::Vector3, 4>& colors);

  std::pair<Id<Handle>, float> nearest_handle(hermite::Vector2 const& position);

  std::pair<Id<HalfEdge>, float> nearest_edge(hermite::Vector2 const& position);

  std::pair<Id<ControlPoint>, float> nearest_point(
      hermite::Vector2 const& position);

  /// Subdivides the selected patch at (t, t2)
  /**
   * Subdivides the patch by three consecutive splits vertical at t and two
   * horizontal splits at t2
   */
  Id<HalfEdge> subdivide(Id<HalfEdge> selected_edge, float t, float t2);

  /// Globally subdivides the mesh at (t, t2)
  /**
   * propagates the splitting of patches along t and t2
   */
  Id<HalfEdge> split_global(Id<HalfEdge> selected_edge, float t, float t2);

  /// Keeps globally splitting neighbouring patches at t relative to top
  /**
   * propagates the splitting of patches along t
   */
  void propagate_split(Id<HalfEdge> top, float t);

  /// Finds value on child that contains t and tranlates to child interval
  /**
   * propagates the splitting of patches along t
   */
  float find_t_on_child(Id<HalfEdge> parent, float t);

  /// finds child that contains t
  /**
   * if no child is found parent is returned
   */
  Id<HalfEdge> find_child_contains_t(Id<HalfEdge> parent, float t);

  /// Find relative t value on child edge
  /**
   *
   */
  float calculate_t_on_parent(Id<HalfEdge> child, float t);

  /// Finds the parent t value of the edge
  /**
   * if its a parent edge next_t t
   * if not it is translated to parent t
   */
  float calculate_next_t(Id<HalfEdge> opp, float t);

  /// returns origin control point
  /**
   */
  Id<ControlPoint> get_origin(Id<HalfEdge> edge);

  /// returs array of the two control points associated with this edge
  /**
   *  This function is implemented in `extrude.cpp`.
   */
  std::array<Id<ControlPoint>, 2> get_control_points(Id<HalfEdge> edge);

  /// extrudes a patch from this boundary edge
  /**
   * will break if not a boundary edge
   * This function is implemented in `extrude.cpp`.
   */
  Id<HalfEdge> extrude(Id<HalfEdge> edge);

  Id<Handle> get_handle_by_number(Id<HalfEdge> edge, int index);

  hermite::Vector2 get_position(Id<ControlPoint> point);

  void move_point(Id<ControlPoint> point, hermite::Vector2 const& delta);

  hermite::Vector2 get_position(Id<Handle> handle);

  void move_handle(Id<Handle> handle, hermite::Vector2 const& delta);

  void set_position(Id<HalfEdge> edge, hermite::Vector2 const& position);

  /// removes all t junctions by propagating their split along the mesh
  /**
   * requires multiple iterations to turn all t edges into full junctions
   * This function is implemented in `splitting.cpp`.
   */
  void remove_t_junctions();

  /// finds t junctions belonging to this edge or if parent is t junction
  /**
   * will return a vector of child t-edges or
   * will return a vector with one item its parent t-edge
   * This function is implemented in `splitting.cpp`.
   */
  std::optional<Id<HalfEdge>> is_parent_t_junction(Id<HalfEdge> edge);

  /// finds t junctions belonging to this edge or if parent is t junction
  /**
   * will return a vector of child t-edges or
   * will return a vector with one item its parent t-edge
   * This function is implemented in `splitting.cpp`.
   */
  std::vector<Id<HalfEdge>> find_t_junctions(Id<HalfEdge> edge);

  /// Splits the patch that top_edge belongs to at point t.
  /**
   * Splits the patch that top_edge belongs to, in the direction
   * parallel to top_edge, at point t.
   * This function is implemented in `splitting.cpp`.
   */
  Id<HalfEdge> split(Id<HalfEdge> top_edge, float t);
  /// Sets a half-edge's color at the origin.
  void set_color(Id<HalfEdge> edge, hermite::Vector3 color);
  void set_color_vertex(Id<HalfEdge> edge, hermite::Vector3 color);
  /// Returns the hermite vector of the given edge.
  hermite::CurveMatrix curve_matrix(Id<HalfEdge> edge) const;
  /**
   * Returns the vector of the  orthogonal curve that would result from
   * splitting edge at point t.
   */
  hermite::CurveMatrix split_curve(Id<HalfEdge> edge, float t) const;
  /**
   * Returns a set of points along the given edge where splitting
   * would create a fully connected junction (snapping).
   * This also considers snap points opposite to the given edge.
   */
  std::set<float> snap_points(Id<HalfEdge> edge) const;

  /**
   * Simplifies the mesh by merging all the neighbouring patches that can be
   * merged without affecting the resulting figure.
   */
  void simplify_mesh(MergeOptions merge_options);
  /**
   * Determines for the two patches neighbouring the specified separator_edge
   * whether they can be merged, based on their continuity.
   */
  bool can_merge_neighbours(Id<HalfEdge> separator_edge,
                            MergeOptions merge_options);
  /**
   * Merges the two patches neighbouring to the given separator_edge when
   * this does not affect the resulting image.
   */
  void merge_neighbours(Id<HalfEdge> separator_edge);

  /// Returns a vector of all of the mesh's patches.
  std::vector<PatchRenderData> patch_data() const;
  /// Returns a vector of all of the mesh's control points.
  std::vector<hermite::Interpolant> point_data() const;
  /// Returns a vector of all of the mesh's tangent handles.
  /**
   * The handles are returned as an array of (start, end),
   * since tangent handles are rendered (also) as line segments.
   */
  std::vector<std::array<hermite::Interpolant, 2>> handle_data() const;
  /// Returns a vector of all of the mesh's boundary curves.
  std::vector<hermite::CurveMatrix> curve_data() const;

  Id<HalfEdge> next(Id<HalfEdge> edge) const;
  Id<HalfEdge> prev(Id<HalfEdge> edge) const;

  bool isChild(Id<HalfEdge> edge) const;
  Interval getInterval(Id<HalfEdge> edge) const;

  void read_from_file(std::string const& file_name);
  void write_to_file(std::string const& file_name) const;

 private:
  /// Creates a new parent half-edge with the given parameters.
  /**
   * It also appropriately links tangent handle backreferences
   * to the returned half-edge.
   */
  Id<HalfEdge> half_edge(Id<ControlPoint> origin,
                         std::array<Id<Handle>, 2> handles,
                         hermite::Vector3 color,
                         hermite::Interpolant twist = {},
                         std::optional<Id<HalfEdge>> twin = std::nullopt);
  /// Creates a new child half-edge with the given parameters.
  /**
   * It also appropriately links tangent handle backreferences
   * to the returned half-edge.
   */
  Id<HalfEdge> half_edge(Id<HalfEdge> parent, Interval interval,
                         std::optional<std::array<Id<Handle>, 2>> handles,
                         hermite::Vector3 color,
                         hermite::Interpolant twist = {},
                         std::optional<Id<HalfEdge>> twin = std::nullopt);

  /// Establishes the neighbor relationship between a and b.
  /**
   * It also copies over the `patch` reference from a to b.
   */
  void connect(Id<HalfEdge> a, Id<HalfEdge> b);
  /// Returns the control matrix of the patch that `edge` belongs to.
  hermite::PatchMatrix patch_matrix(Id<HalfEdge> edge) const;
  /**
   * Returns the four interpolants of a half-edge,
   * in the order `[m0, m0uv, m0v, m1v]`.
   */
  std::array<hermite::Interpolant, 4> edge_tangents(Id<HalfEdge> edge) const;
  /// Compute the boundary information of the given edge.
  /**
   * Boundary information includes the parent curve
   * and the segments to be used with tessellation.
   */
  EdgeBoundary boundary(Id<HalfEdge> edge) const;

  /// Compute the boundary information of the given patch.
  PatchBoundary boundary(Patch const& patch) const;

  /// Returns the two endpoints of a handle line.
  std::array<hermite::Interpolant, 2> endpoints(const Handle& handle) const;

  //=========================================================================
  // The following functions are related to splitting,
  // and thus implemented in splitting.cpp.
  //=========================================================================

  /// Creates a child edge parented to the given edge.
  /**
   * The edge will be at parameter value 0 along the parent.
   * If the edge is already a child, this function will simply adjust
   * the edge's length, without creating a new one.
   * @param edge: the edge to create a child from.
   * @param length: the length of the child relative to the parent.
   */
  Id<HalfEdge> child_from(Id<HalfEdge> edge, float length);
  /// Connects the given half-edge to the control point `origin`.
  /**
   * This process transforms a child edge to a parent edge.
   */
  void to_parent(Id<HalfEdge> edge, Id<ControlPoint> origin);

  /**
   * Returns the  half-edges originating from splitting the given edge,
   * in the order [left, middle, right].
   */
  std::array<Id<HalfEdge>, 3> junction(Id<HalfEdge> edge, float t,
                                       hermite::Interpolant twist,
                                       std::array<Id<Handle>, 2> ortho_handles);
  /**
   * Creates a "connected" junction (where the mid half-edge is connected to a
   * control point) from the given parameters.
   */
  std::array<Id<HalfEdge>, 3> connected_junction(
      Id<HalfEdge> edge, float t, hermite::Interpolant twist,
      std::array<Id<Handle>, 2> ortho_handles);

  /**
   * Creates a T-junction (where the mid half-edge lays along a parent curve)
   * at parameter value t along the given edge.
   */
  std::array<Id<HalfEdge>, 3> t_junction(
      Id<HalfEdge> edge, float t, hermite::Interpolant twist,
      std::array<Id<Handle>, 2> ortho_handles);

  /**
   * If the given edge is sufficiently close to its twin, this
   * function will snap them together in a connected junction.
   */
  void snap_to_twin(Id<HalfEdge> edge);

  /**
   * Creates a "connected" junction (where the mid half-edge is connected to a
   * control point) from the given parameters. This is to be used when snapping,
   * compared to `connected_junction.`
   */
  std::array<Id<HalfEdge>, 2> parent_junction(
      Id<HalfEdge> edge, float t, Id<ControlPoint> mid_point,
      hermite::Vector3 color, hermite::Interpolant twist,
      std::array<Id<Handle>, 2> parallel_handles);

  /**
   * Reparents all edges in the range [begin, end) to `parent`, and adjusts
   * their parameters according to `new_interval`.
   * When snapping two twin edges, moves all of their siblings to their
   * appropriate points and lengths along their new parents.
   */
  void reparametrize(ChildIterator begin, ChildIterator end,
                     Id<HalfEdge> parent, Interval new_interval);

  /// Multiplies the tangent handles connected to the given edge by t.
  void scale_tangents(Id<HalfEdge> edge, float t);

  /// Multiplies the twist vectors connected to the given edge by t.
  void scale_twists(Id<HalfEdge> edge, float t);

  //=========================================================================
  // The following functions are related to merging,
  // and thus implemented in merge-conditions.cpp and merging.cpp.
  //=========================================================================

  /**
   * Returns a potential splitting factor `t` based on the specified left and
   * right matrices.
   */
  float potential_splitting_factor(hermite::PatchMatrix left,
                                   hermite::PatchMatrix right);

  /**
   * Returns the relative (parent) interval lengths using the splitting factor
   * t, namely:
   *  - The length of the left interval, relative to the adjacent interval
   *  - The length of the adjacent interval, relative to the left interval
   * This can be used to find the interval necessary to reparametrize children
   * of the left and adjacent intervals to their merged parent interval.
   * See also merge_parents(...).
   */
  std::array<float, 2> relative_interval_lengths(Id<HalfEdge> left,
                                                 Id<HalfEdge> adjacent,
                                                 float t);

  /**
   * Calculates the control matrix of the patch acquired when merging left
   * and right. It assumes they can be merged; if not, nonsense is returned.
   */
  hermite::PatchMatrix parent(hermite::PatchMatrix left,
                              hermite::PatchMatrix right, float t);

  /**
   * Merges the outside of left and adjacent when necessary. Returns the
   * twin-edge when one exists.
   */
  std::optional<Id<HalfEdge>> merge_outside(Id<HalfEdge> left,
                                            Id<HalfEdge> adjacent,
                                            float relative_left_length,
                                            float relative_adjacent_length);

  /**
   * Merges the inside of the specified left and adjacent edges.
   * If the junction is a connected-junction, it first transforms it to a
   * t-junction. Then the left- and adjacent edges are merged.
   * Returns references to the left child, left_parent and, for removal,
   * the orthogonal edge and origin.
   */
  std::tuple<Id<HalfEdge>, Id<HalfEdge>, Id<HalfEdge>,
             std::optional<Id<ControlPoint>>>
  merge_inside(Id<HalfEdge> left, Id<HalfEdge> adjacent,
               hermite::CurveMatrix edge_row, float relative_left_length,
               float relative_adjacent_length);

  /**
   * Returns a parent-child pair using the following logic:
   * - supplied edge is parentable: create child and return it with parent.
   * - supplied edge is not parentable: return it with it's parent.
   * This assumes the supplied edge does not have children yet.
   */
  std::array<Id<HalfEdge>, 2> parent_child(Id<HalfEdge> edge);

  /**
   * Merge left_parent & adjacent_parent by appending children of
   * adjacent_parent to left_parent and merging parent edges.
   * Returns references to the left_parent and the ControlPoint of the
   * orthogonal edge, the latter so that it can be removed later.
   */
  std::pair<Id<HalfEdge>, Id<ControlPoint>> merge_parents(
      Id<HalfEdge> left_parent, Id<HalfEdge> adjacent_parent,
      Interval new_left_interval, Interval new_adjacent_interval);

  /**
   * Removes an edge. The handles of the edge are removed when it does not have
   * a twin.
   */
  void remove_edge(Id<HalfEdge> edge);

  /**
   * Merges the specified left- and adjacent edges. It assumes the junction
   * between the two edges is a t-junction.
   */
  Id<HalfEdge> merge_adjacent_siblings(Id<HalfEdge> left,
                                       Id<HalfEdge> adjacent);

  /**
   * Simplifies a parent-child pair by merging them, if possible.
   * You can supply either a child or parent. It returns the edge
   * itself when simplification was not possible, or the parent
   * after simplification.
   */
  Id<HalfEdge> merge_parent_with_child(Id<HalfEdge> parent_or_child);

  /**
   * Removes a parent edge and all its first descendents, including the handles
   * (i.e. first child, first child of first child, etc.).
   */
  void remove_edge_parent_child(Id<HalfEdge> parent);

  /**
   * Removes the handles of an edge.
   */
  void remove_handles(Id<HalfEdge> edge);

  /**
   * Sets the twin of the specified parent and all its children to the
   * specified twin.
   */
  void set_twin(Id<HalfEdge> parent, Id<HalfEdge> twin);

  /**
   * Update the handle references of twin to the handles of source,
   * removing the old handles of twin.
   */
  void share_handles(Id<HalfEdge> source, Id<HalfEdge> twin);

  /**
   * Update the twists of edge and it's neighbour to the specified twists.
   */
  void update_twists(Id<HalfEdge> edge,
                     std::array<hermite::Interpolant, 2> new_twists);

  void read_points(std::istream& input, int num_points);
  void read_handles(std::istream& input, int num_handles);
  void read_patches(std::istream& input, int num_patches);
  void read_edges(std::istream& input, int num_edges);

  void write_points(std::ostream& output) const;
  void write_handles(std::ostream& output) const;
  void write_patches(std::ostream& output) const;
  void write_edges(std::ostream& output) const;

  std::vector<float> rgb_data() const;
  std::vector<float> rgbxy_data() const;

  std::vector<hermite::Vector3> get_palette() const;
  std::vector<hermite::Vector3> get_recolored(
      std::vector<hermite::Vector3> palette) const; // TODO: get_weights.

  void recolor(std::vector<hermite::Vector3> rgb);

  Storage<Handle> handles;
  Storage<HalfEdge> edges;
  Storage<Patch> patches;
  Storage<ControlPoint> points;
};

} // namespace gmt
