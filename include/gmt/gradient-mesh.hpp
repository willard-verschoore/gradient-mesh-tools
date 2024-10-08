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
  Id<HalfEdge> find_child_contains_t(Id<HalfEdge> parent, float t) const;

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

  /// Determines the continuity across the two patches neighboring \c separator.
  /**
   * Continuity up to C4 is checked. The returned integer indicates the degree
   * of continuity. If the patches are discontinuous or \c separator is invalid,
   * then -1 is returned. The separator edge is invalid if it is a boundary edge
   * and therefore does not neighbor two patches. A parent edge with children is
   * also invalid as this means that at least one of the neighboring patches is
   * split.
   *
   * @param separator: The edge separating the two patches for which the
   * continuity will be found.
   * @param options: The tolerances to use when checking for continuity.
   * @return An integer from -1 to 4 indicating the degree of continuity.
   */
  int patch_continuity(Id<HalfEdge> separator, MergeOptions options) const;

  /// Determines the continuity across every pair of neighboring patches.
  /**
   * Continuity up to C4 is checked. The returned integers indicate the degrees
   * of continuity. If two patches are discontinuous the corresponding degree is
   * -1.
   *
   * @param options: The tolerances to use when checking for continuity.
   * @return An array of integers from -1 to 4 indicating the degrees of
   * continuity.
   */
  std::vector<int> patch_continuities(MergeOptions options) const;

  /// Removes children from parents which have only one child.
  void remove_only_children();

  /// Returns the color at the leaf of the tree below \c edge where t = 0.
  hermite::Vector3 get_origin_color(Id<HalfEdge> edge) const;

  /// Sets the color for every node of the tree below \c edge where t = 0.
  void set_origin_color(Id<HalfEdge> edge, hermite::Vector3 color);

  /// Ensures that the color of an inactive edge is the active color at t = 0.
  void fix_inactive_origin_colors();

  /// Returns a count of the number of handles in the mesh.
  int handle_count() const;

  /// Returns a count of the number of edges in the mesh.
  int edge_count() const;

  /// Returns a count of the number of patches in the mesh.
  int patch_count() const;

  /// Returns a count of the number of control points in the mesh.
  int point_count() const;

  std::vector<PatchRenderData> render_data() const;
  /// Returns a vector of all of the mesh's patches.
  std::vector<hermite::PatchMatrix> patch_data() const;
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

  /// Reads the data from the given PatchMatrix objects into the mesh.
  /**
   * Ideally, after calling this function a call to patch_data() would return an
   * array identical to \c data. There are some restrictions to this however.
   *
   * This is because an array of patch matrices contains redundant data. For
   * example, the coordinates of a control point will occur in the matrix of
   * every adjacent patch. If \c data contains contradictory coordinates for
   * such a control point only the values in the last patch matrix will actually
   * be used.
   *
   * A similar situation occurs for patches resulting from splits where at least
   * one side is a child edge. Some properties of such an edge are determined by
   * its parent. If \c data contains new values for these properties they will
   * be ignored. Changes to the parent edge will be applied and propagated to
   * its children.
   *
   * @param data The PatchMatrix objects to read from.
   */
  void read_patch_data(std::vector<hermite::PatchMatrix> const& data,
                       bool create_tangents);

  void read_from_file(std::string const& file_name);
  void write_to_file(std::string const& file_name) const;

  static const int MAX_SAMPLING_DENSITY;

  /// Samples colors from each patch in the mesh.
  /**
   * Makes use of a patch's Hermite function H(u, v) to sample colors from each
   * patch. The function samples using uniformly spaced coordinates in uv-space
   * (i.e. from [0, 1]^2). The \c density parameter determines the number of
   * samples to take in one axis. A density of 1 corresponds to sampling only
   * the patch corners for a total of 4 colors per patch. In general (density +
   * 1)^2 colors will be sampled per patch.
   *
   * @param density The number of samples to take per patch in each axis.
   * @return An array of colors sampled from the mesh's patches.
   */
  std::vector<hermite::Vector3> sample_colors(int density = 1) const;

  /// Extracts the mesh's control point colors as a flat array of floats.
  /**
   * The colors are taken from the mesh's patch matrices, either in Hermite or
   * Bezier form, as specified by the \c bezier parameter.
   *
   * The Hermite form is used by the library internally and the form provided
   * through the patch_data() function. Only the 4 corners of a Hermite matrix
   * are proper RGB(XY) control points though, the other elements are vectors in
   * RGB(XY) space. This function thus only returns 4 colors per patch in the
   * Hermite case.
   *
   * For Bezier patch matrices each element corresponds to a control point and
   * is thus a valid RGB(XY) position. When \c bezier is set this function
   * therefore returns 16 colors per patch.
   *
   * The \c inactive parameter specifies whether data from the mesh's inactive
   * edges is also included. An inactive edge is an edge which has child edges
   * running in parallel with it. An inactive edge is not directly adjacent to
   * any patch, instead it controls its parallel children, which are each
   * adjacent to their own patch. This means that an inactive edge's data is not
   * included when retrieving patch data but it may still be relevant due to the
   * effect it has on its children.
   *
   * The data of a curve corresponding to an inactive edge consists of two
   * endpoints and two tangents. Similarly to the data from the patch matrices,
   * the tangents are only valid control points when the curve is converted to
   * Bezier form. Therefore, when \c bezier is set each inactive edge adds 4
   * colors and when it is not set only the 2 endpoint colors are included.
   *
   * @param bezier Whether to get the RGB data from the Bezier form of the
   * mesh's patches.
   * @param inactive Whether to get the RGB data from the mesh's inactive edges.
   * @return An array of floats where each group of three consecutive elements
   * represents an RGB color of a control point in the mesh.
   */
  std::vector<float> rgb_data(bool bezier, bool inactive) const;

  /// Extracts the mesh's control point RGBXY data as a flat array of floats.
  /**
   * The RGBXY positions are taken from the mesh's patch matrices, either in
   * Hermite or Bezier form, as specified by the \c bezier parameter.
   *
   * The Hermite form is used by the library internally and the form provided
   * through the patch_data() function. Only the 4 corners of a Hermite matrix
   * are proper RGBXY control points though, the other elements are vectors in
   * RGBXY space. This function thus only returns 4 RGBXY positions per patch in
   * the Hermite case.
   *
   * For Bezier patch matrices each element corresponds to a control point and
   * is thus a valid RGBXY position. When \c bezier is set this function
   * therefore returns 16 RGBXY positions per patch.
   *
   * The \c inactive parameter specifies whether data from the mesh's inactive
   * edges is also included. An inactive edge is an edge which has child edges
   * running in parallel with it. An inactive edge is not directly adjacent to
   * any patch, instead it controls its parallel children, which are each
   * adjacent to their own patch. This means that an inactive edge's data is not
   * included when retrieving patch data but it may still be relevant due to the
   * effect it has on its children.
   *
   * The data of a curve corresponding to an inactive edge consists of two
   * endpoints and two tangents. Similarly to the data from the patch matrices,
   * the tangents are only valid control points when the curve is converted to
   * Bezier form. Therefore, when \c bezier is set each inactive edge adds 4
   * RGBXY positions and when it is not set only the 2 endpoint RGBXY positions
   * are included.
   *
   * @param bezier Whether to get the RGBXY data from the Bezier form of the
   * mesh's patches.
   * @param inactive Whether to get the RGBXY data from the mesh's inactive
   * edges.
   * @return An array of floats where each group of five consecutive elements
   * represents an RGB color and XY position of a control point in the mesh.
   */
  std::vector<float> rgbxy_data(bool bezier, bool inactive) const;

  /// Extracts a palette for the gradient mesh using the RGB convex hull.
  /**
   * The palette is found by taking the RGB-space convex hull of colors
   * sampled from the mesh. The number of colors to sample is determined by \c
   * sampling_density. See the documentation of sample_colors() for more
   * details. The vertices of the convex hull are the palette colors. The convex
   * hull is simplified to the number of vertices specified by \c target_size.
   *
   * Besides the palette color, this function also returns the palette indices.
   * The indices come in pairs where each pair indicates a connection between
   * two palette colors.
   *
   * @param target_size The target number of palette colors.
   * @param sampling_density The density with which to sample colors from the
   * mesh.
   * @return A pair of palette colors and indices.
   */
  std::pair<std::vector<hermite::Vector3>, std::vector<uint32_t>> get_palette(
      size_t target_size, int sampling_density = 1) const;

  /// Simplifies a palette by attempting to remove one color.
  /**
   * Tries to remove a vertex from the palette hull by collapsing an edge into a
   * single point. The edge collapse that leads to the smallest added volume is
   * selected. There may be no possible edge collapses, in which case this
   * function returns the original input palette.
   *
   * Besides the new palette colors, this function also returns new palette
   * indices. The indices come in pairs where each pair indicates a connection
   * between two palette colors.
   *
   * TODO: Consider making this function static or moving it out of the
   * GradientMesh class all together. It does not use any mesh data.
   *
   * @param palette The palette to simplify. Generally this should be the result
   * of get_palette().
   * @param sampling_density The density with which to sample colors from the
   * mesh.
   * @return A simplified pair of palette colors and indices.
   */
  std::pair<std::vector<hermite::Vector3>, std::vector<uint32_t>>
  simplify_palette(std::vector<hermite::Vector3> const& palette) const;

  /// Optimizes a palette's vertex positions.
  /**
   * Tries to find a trade-off between minimizing reconstruction and
   * representation loss. Reconstruction loss is caused by mesh colors outside
   * of the palette hull. Representation loss is caused by palette colors far
   * away from any mesh colors. For more details on this method, see the
   * original paper: <a
   * href=https://www.doi.org/10.1111/cgf.13812>doi.org/10.1111/cgf.13812</a>.
   *
   * To obtain more detailed results for the reconstruction and representation
   * loss one can increase the \c sampling_density parameter. See the
   * documentation of sample_colors() for more details. Note that this may
   * significantly slow down the method.
   *
   * @param palette The palette to optimize. Generally this should be the result
   * of get_palette().
   * @param lambda The lambda parameter in the energy function from the original
   * paper. Reconstruction loss is weighed more heavily for higher values.
   * @param sampling_density The density with which to sample colors from the
   * mesh.
   * @return An optimized version of the input palette.
   */
  std::vector<hermite::Vector3> optimize_palette(
      std::vector<hermite::Vector3> const& palette, double lambda = 10.0f,
      int sampling_density = 1) const;

  /// The different types of palette weights that are supported.
  enum class WeightType
  {
    /**
     * Weights determined based on the RGBXY-space geometry the mesh. For more
     * details on this method, see the original paper: <a
     * href=https://doi.org/10.1145/3272127.3275054>doi.org/10.1145/3272127.3275054</a>.
     */
    RGBXY,

    /**
     * Mean Value Coordinates computed from the triangular mesh corresponding to
     * a palette's convex hull. See the following paper for more details: <a
     * href=https://www.doi.org/10.1111/cgf.13812>doi.org/10.1111/cgf.13812</a>.
     */
    MVC
  };

  /// Finds weights for the palette colors which reproduce each mesh color.
  /**
   * Finds weight vectors that specify linear combinations of the palette colors
   * which reproduce the colors in the mesh.
   *
   * @param palette The palette colors to determine the weights for. Generally
   * this should be the result of get_palette().
   * @param type The type of weights to get. See the documentation of
   * GradientMesh::WeightType for more details.
   * @param bezier Whether to get the mesh's RGBXY positions using the Bezier
   * form of the mesh's patch matrices. See the documentation of rgbxy_data()
   * for more details.
   * @param inactive Whether to include inactive edges when retrieving the
   * mesh's RGBXY positions. See the documentation of rgbxy_data() for more
   * details.
   * @return weights The weights for the palette colors which reproduce each
   * mesh color. Every P consecutive elements specify the weights for one point
   * in the mesh, where P is the number of palette colors.
   */
  std::vector<float> get_weights(std::vector<hermite::Vector3> const& palette,
                                 WeightType type, bool bezier,
                                 bool inactive) const;

  /// Applies weights to a palette to obtain a recolored version of the mesh.
  /**
   * @param weights Weights that specify the contribution of the palette colors
   * for each point in the mesh. Generally this should be the result of
   * get_weights().
   * @param palette The palette with which to recolor the mesh. If this is the
   * result of get_palette() the colors of the mesh do not change.
   * @param bezier Whether the weights were obtained using the Bezier form of
   * the mesh's patch matrices. Should be equivalent to the value passed to the
   * matching get_weights() call.
   * @param inactive Whether the weights were obtained using data from the
   * mesh's inactive edges. Should be equivalent to the value passed to the
   * matching get_weights() call.
   * @param create_tangents Whether to create new tangent handles when
   * recoloring edges that do not have any. This gives better results but may
   * lead to issues when splitting later on. Should only be enabled when using
   * Bezier weights.
   */
  void recolor(std::vector<float> const& weights,
               std::vector<hermite::Vector3> const& palette, bool bezier,
               bool inactive, bool create_tangents = false);

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
                         bool recolored, hermite::Vector3 color,
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
                                   hermite::PatchMatrix right) const;

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

  /// Reads the given data into \c edge.
  /**
   * Here \c data should be of the form [origin, twist, handle0, handle1]. The
   * color component of the provided origin and the entire twist vector are
   * always read into the edge. If the edge is a parent the origin coordinates
   * and handle vectors are also read.
   *
   * If \c edge is a child the origin coordinates are ignored. The handle
   * vectors are only read if the child edge has handles already, unless the \c
   * create_tangents option is specified. In that case new handles will be
   * created.
   *
   * @param edge The edge to read the data into.
   * @param data The Interpolant objects to read into the edge.
   * @param create_tangents Whether to create new tangent handles for child
   * edges that do not have them.
   */
  void read_edge(Id<HalfEdge> edge,
                 std::array<hermite::Interpolant, 4> const& data,
                 bool create_tangents);

  /// Reads the data from the given CurveMatrix object into \c edge.
  /**
   * Ideally, after calling this function a call to curve_matrix() for \c edge
   * would return a CurveMatrix identical to \c data. There are some limitations
   * to this however.
   *
   * This is because child edges determine their origin coordinates based on
   * their parent. A child may also have tangent handles determined by their
   * parent. For the latter it is possible to forcefully create new tangents to
   * match those in \c data by setting the \c create_tangents parameter.
   *
   * Note that this function is similar to, but slightly different from, the
   * read_edge() function. This function reads in a CurveMatrix containing the
   * two endpoints of a curve and its two tangent handles. Meanwhile,
   * read_edge() expects only the endpoint directly stored in \c edge, the same
   * two tangent handles, and a twist vector. Broadly speaking, this function
   * matches the mathematical definition of a Hermite curve while read_edge()
   * matches the way a HalfEdge is defined in the code.
   *
   * @param edge The edge to read the data into.
   * @param data The CurveMatrix object to read from.
   * @param create_tangents Whether to create new tangent handles for child
   * edges that do not have them.
   */
  void read_curve_matrix(Id<HalfEdge> edge, hermite::CurveMatrix const& matrix,
                         bool create_tangents);

  /// Reads the data from the given PatchMatrix object into \c patch.
  /**
   * Ideally, after calling this function a call to patch_matrix() for \c patch
   * would return a PatchMatrix identical to \c data. There are some limitations
   * to this however.
   *
   * This is because child edges in a patch determine their origin coordinates
   * based on their parent. A child may also have tangent handles determined by
   * their parent. For the latter it is possible to forcefully create new
   * tangents to match those in \c data by setting the \c create_tangents
   * parameter.
   *
   * @param patch The patch to read the data into.
   * @param data The PatchMatrix object to read from.
   * @param create_tangents Whether to create new tangent handles for child
   * edges that do not have them.
   */
  void read_patch_matrix(Patch const& patch, hermite::PatchMatrix const& matrix,
                         bool create_tangents);

  /// Finds the color on the other side of the edge at the origin of \c edge.
  hermite::Vector3 find_color_opposite_origin(Id<HalfEdge> edge) const;

  /// Finds t-junctions whose color matches the other side of the edge.
  std::vector<Id<HalfEdge>> find_continuous_t_junctions() const;

  /// Ensures that each junction's color matches the other side of the edge.
  void fix_continuous_t_junctions(std::vector<Id<HalfEdge>> const& junctions);

  void read_points(std::istream& input, int num_points);
  void read_handles(std::istream& input, int num_handles);
  void read_patches(std::istream& input, int num_patches);
  void read_edges(std::istream& input, int num_edges);

  void write_points(std::ostream& output) const;
  void write_handles(std::ostream& output) const;
  void write_patches(std::ostream& output) const;
  void write_edges(std::ostream& output) const;

  size_t get_index(Id<Handle> handle) const;
  size_t get_index(Id<HalfEdge> handle) const;
  size_t get_index(Id<Patch> handle) const;
  size_t get_index(Id<ControlPoint> handle) const;

  Storage<Handle> handles;
  Storage<HalfEdge> edges;
  Storage<Patch> patches;
  Storage<ControlPoint> points;
};

} // namespace gmt
