#include <cmath>

#include "gradient-mesh.hpp"

/**
 * Calculates the original patch control matrix in hermite form,
 * based on the left patch and the splitting value `t`.
 */
inline hermite::PatchMatrix parent_left(hermite::PatchMatrix left, float t)
{
  auto mat = hermite::matrix<4, 4>({1.0f, 0.0f, 0.0f, 0.0f,             //
                                    0.0f, 1.0f / t, 0.0f, 0.0f,         //
                                    0.0f, 0.0f, 1.0f / pow(t, 2), 0.0f, //
                                    0.0f, 0.0f, 0.0f, 1.0f / pow(t, 3)});
  auto s_Linv = hermite::n_inv * mat * hermite::n;
  return left * s_Linv.transposed();
}

/**
 * Calculates the control matrix of the original patch in hermite form,
 * based on the right patch and the splitting value `t`.
 */
inline hermite::PatchMatrix parent_right(hermite::PatchMatrix right, float t)
{
  auto t2 = t * t;
  auto t3 = t2 * t;
  auto tf = 1.0f - t;
  auto tf2 = pow(tf, 2);
  auto tf3 = pow(tf, 3);
  auto mat =
      hermite::matrix<4, 4>({1.0f, -t / tf, t2 / tf2, -t3 / tf3,          //
                             0.0f, 1.0f / tf, -2 * t / tf2, 3 * t2 / tf3, //
                             0.0f, 0.0f, 1.0f / tf2, -3 * t / tf3,        //
                             0.0f, 0.0f, 0.0f, 1.0f / tf3});
  auto s_Rinv = hermite::n_inv * mat * hermite::n;
  return right * s_Rinv.transposed();
}

/**
 * Return the specified CurveMatrix, but flipped.
 */
inline hermite::CurveMatrix flip(hermite::CurveMatrix mat)
{
  return hermite::matrix<4, 1>({mat(0, 3), -mat(0, 2), -mat(0, 1), mat(0, 0)});
}

// Simple implementation, multiple iterations are necessary.
void GradientMesh::simplify_mesh(MergeOptions merge_options)
{
  int edge_counter = 0;
  int merge_counter = 0;
  auto it = patches.begin();
  while (it != patches.end())
  {
    auto edge = it->side;
    for (int i = 0; i < 4; ++i)
    {
      if (!edges.contains(edge))
      {
        break;
      }
      auto next = edges[edge].next;
      edge_counter++;

      if (can_merge_neighbours(edge, merge_options))
      {
        merge_neighbours(edge);
        it = patches.begin();
        merge_counter++;
      }

      edge = next;
    }
    ++it;
  }
}

void GradientMesh::merge_neighbours(Id<HalfEdge> separator_edge)
{
  auto middle_left = separator_edge;
  auto bottom_left = edges[middle_left].next;
  auto left = edges[bottom_left].next;
  auto top_left = edges[left].next;

  auto middle_right = edges[separator_edge].twin.value();
  auto top_right = edges[middle_right].next;

  // Since the top_right is determined using adjacent_next, it can be a parent.
  // In that case, we should use the child, otherwise determining the patch
  // matrix is very difficult.
  top_right = first_child(top_right, edges);

  auto right = edges[top_right].next;
  auto bottom_right = edges[right].next;

  auto patch_left = edges[left].patch;
  auto patch_right = edges[right].patch;

  // Calculate potential splitting factor
  auto matrix_left = patch_matrix(top_left);
  auto matrix_right = patch_matrix(top_right);
  auto t = potential_splitting_factor(matrix_left, matrix_right);

  // Calculate hermite control matrix of new patch.
  auto hermite_parent = parent(matrix_left, matrix_right, t);

  // Calculate new relative intervals
  auto [relative_top_left_length, relative_top_adjacent_length] =
      relative_interval_lengths(top_left, top_right, t);
  auto [relative_bottom_left_length, relative_bottom_adjacent_length] =
      relative_interval_lengths(bottom_right, bottom_left, 1.0f - t);

  // Handle outside edges
  auto top_twin = merge_outside(top_left, top_right, relative_top_left_length,
                                relative_top_adjacent_length);
  auto bottom_twin =
      merge_outside(bottom_right, bottom_left, relative_bottom_left_length,
                    relative_bottom_adjacent_length);

  // Handle inside edges
  auto [top_left_child, top_left_parent, top_orthogonal,
        top_orthogonal_origin] =
      merge_inside(top_left, top_right, hermite::row(hermite_parent, 0),
                   relative_top_left_length, relative_top_adjacent_length);
  auto [bottom_right_child, bottom_right_parent, bottom_orthogonal,
        bottom_orthogonal_origin] =
      merge_inside(
          bottom_right, bottom_left, flip(hermite::row(hermite_parent, 3)),
          relative_bottom_left_length, relative_bottom_adjacent_length);

  // Simplify edges
  top_left = merge_parent_with_child(top_left_child);
  bottom_right = merge_parent_with_child(bottom_right_child);

  // Remove orthogonal edges & points.
  remove_edge_parent_child(top_orthogonal);
  if (top_orthogonal_origin.has_value())
  {
    points.remove(top_orthogonal_origin.value());
  }
  remove_edge_parent_child(bottom_orthogonal);
  if (bottom_orthogonal_origin.has_value())
  {
    points.remove(bottom_orthogonal_origin.value());
  }

  // Simplify twin edges and repair twins & handles.
  if (top_twin.has_value())
  {
    auto actual_top_twin = merge_parent_with_child(top_twin.value());
    set_twin(top_left_parent, actual_top_twin);
    set_twin(actual_top_twin, top_left_parent);

    share_handles(top_left_parent, actual_top_twin);
  }
  if (bottom_twin.has_value())
  {
    auto actual_bottom_twin = merge_parent_with_child(bottom_twin.value());
    set_twin(bottom_right_parent, actual_bottom_twin);
    set_twin(actual_bottom_twin, bottom_right_parent);

    share_handles(bottom_right_parent, actual_bottom_twin);
  }

  // Connect border edges of new patch.
  ::connect(left, patch_left, edges, patches);
  connect(left, top_left);
  connect(top_left, right);
  connect(right, bottom_right);
  connect(bottom_right, left);

  // Update twist vectors.
  auto m0uv = hermite_parent(1, 1);
  auto m1uv = -hermite_parent(1, 2);
  auto m2uv = -hermite_parent(2, 1);
  auto m3uv = hermite_parent(2, 2);
  update_twists(top_left, {m0uv, m1uv});
  update_twists(bottom_right, {m3uv, m2uv});

  patches.remove(patch_right);
}

hermite::PatchMatrix GradientMesh::parent(hermite::PatchMatrix left,
                                          hermite::PatchMatrix right, float t)
{
  auto mat = parent_left(left, t) + parent_right(right, t);
  mat *= 0.5;
  return mat;
}

std::array<float, 2> GradientMesh::relative_interval_lengths(
    Id<HalfEdge> left, Id<HalfEdge> adjacent, float t)
{
  auto left_length = length(edges[left].relative_interval());
  auto adjacent_length = length(edges[adjacent].relative_interval());
  auto relative_left_length =
      (adjacent_length / (1.0f - t) - adjacent_length) / left_length;
  auto relative_adjacent_length =
      (left_length / t - left_length) / adjacent_length;
  return {relative_left_length, relative_adjacent_length};
}

std::optional<Id<HalfEdge>> GradientMesh::merge_outside(
    Id<HalfEdge> left, Id<HalfEdge> adjacent, float relative_left_length,
    float relative_adjacent_length)
{
  // Outer edges, nothing to merge.
  if (!edges[left].twin.has_value() || !edges[adjacent].twin.has_value())
  {
    return std::nullopt;
  }

  auto left_twin = edges[left].twin.value();
  auto adjacent_twin = edges[adjacent].twin.value();

  // Other side is t-junction or has no junction.
  if (left_twin == adjacent_twin)
  {
    return std::make_optional(adjacent_twin);
  }

  // Ensure both twins are a parent-child combo.
  if (left_twin == last_child(left_twin, edges))
  {
    parent_child(left_twin);
  }
  if (adjacent_twin == first_child(adjacent_twin, edges))
  {
    parent_child(adjacent_twin);
  }

  // Merge parent-child combinations.
  auto [left_parent, unused] = merge_parents(
      adjacent_twin, left_twin, Interval{0.0f, 1.0f + relative_left_length},
      Interval{0.0f - relative_adjacent_length, 1.0f});
  return std::make_optional(left_parent);
}

std::tuple<Id<HalfEdge>, Id<HalfEdge>, Id<HalfEdge>,
           std::optional<Id<ControlPoint>>>
GradientMesh::merge_inside(Id<HalfEdge> left, Id<HalfEdge> adjacent,
                           hermite::CurveMatrix edge_row,
                           float relative_left_length,
                           float relative_adjacent_length)
{
  auto [left_parent, left_child] = parent_child(left);
  auto [adjacent_parent, adjacent_child] = parent_child(adjacent);
  std::optional<Id<ControlPoint>> orthogonal_origin = std::nullopt;

  // Junction is connected junction, transform to t-junction
  if (left_parent != adjacent_parent)
  {
    auto [new_left_parent, new_orthogonal_origin] =
        merge_parents(left_parent, adjacent_parent,
                      Interval{0.0f, 1.0f + relative_adjacent_length},
                      Interval{0.0f - relative_left_length, 1.0f});
    left_parent = new_left_parent;
    orthogonal_origin = std::make_optional(new_orthogonal_origin);

    // Update handles & twists.
    auto new_patch_start = edges[left_child].relative_interval().start;
    auto new_patch_end = edges[adjacent_child].relative_interval().end;
    auto new_patch_length = new_patch_end - new_patch_start;
    auto left_handle = hermite::parallel_derivative(
                           edge_row, -new_patch_start / new_patch_length) /
                       new_patch_length;
    auto adjacent_handle =
        hermite::parallel_derivative(
            edge_row, 1.0f + (1.0f - new_patch_end) / new_patch_length) /
        new_patch_length;
    auto next = edges[left_parent].next;
    auto h = edges[left_parent].handles();
    assert(h.has_value());
    auto [start, end] = h.value();
    handles[start].tangent = left_handle;
    handles[end].tangent = -adjacent_handle;
    handles[end].edge = next;

    edges[left_parent].twist /= new_patch_length;
    edges[next].twist /= new_patch_length;
  }

  auto orthogonal_edge = merge_adjacent_siblings(left_child, adjacent_child);

  return std::make_tuple(left_child, left_parent, orthogonal_edge,
                         orthogonal_origin);
}

std::array<Id<HalfEdge>, 2> GradientMesh::parent_child(Id<HalfEdge> edge)
{
  if (edges[edge].is_parentable())
  {
    auto child = child_from(edge, 1.0f);

    // Connect previous edge(s)
    auto prev = edges[edge].prev;
    auto prev_parent = eligible_parent(prev, edges);
    // Only connect prev_parent when prev is the last child.
    if (prev_parent != prev && edges[prev_parent].next == edge)
    {
      connect(prev_parent, child);
    }
    connect(prev, child);

    // Connect next edge and patch.
    connect(child, edges[edge].next);
    ::connect(child, edges[edge].patch, edges, patches);
    return {edge, child};
  }
  return {eligible_parent(edge, edges), edge};
}

std::pair<Id<HalfEdge>, Id<ControlPoint>> GradientMesh::merge_parents(
    Id<HalfEdge> left_parent, Id<HalfEdge> adjacent_parent,
    Interval new_left_interval, Interval new_adjacent_interval)
{
  auto left_children = children(left_parent, edges);
  auto adjacent_children = children(adjacent_parent, edges);

  // Acquire necessary references.
  auto last_left_child = std::prev(left_children.end())->second;
  auto first_adjacent_child = adjacent_children.begin()->second;

  // Append orthogonal edge to left children.
  auto orthogonal_edge = eligible_parent(edges[last_left_child].next, edges);
  auto orthogonal_origin = std::get<Parent>(edges[orthogonal_edge].kind).origin;
  auto last_left_child_end = edges[last_left_child].interval().end;
  edges[orthogonal_edge].kind =
      Child{left_parent, Interval{last_left_child_end, last_left_child_end},
            edges[orthogonal_edge].handles()};

  // Adopt left children, fixing the intervals.
  reparametrize(left_children.begin(), left_children.end(), left_parent,
                new_left_interval);

  // Ensure adjacent parent is no longer referenced.
  auto orthogonal_twin = edges[orthogonal_edge].twin;
  auto prev = edges[first_adjacent_child].prev;
  assert(orthogonal_twin.has_value());
  assert(eligible_parent(prev, edges) == orthogonal_twin.value());
  connect(orthogonal_twin.value(), first_adjacent_child);
  connect(prev, first_adjacent_child);

  // Adopt adjacent children, fixing the intervals.
  reparametrize(adjacent_children.begin(), adjacent_children.end(), left_parent,
                new_adjacent_interval);

  // Restore original leftmost child & link left_parent.
  edges[left_parent].leftmost_child = left_children.begin()->second;
  edges[left_parent].next = edges[adjacent_parent].next;

  // Remove no longer used adjacent parent.
  remove_edge(adjacent_parent);

  return std::make_pair(left_parent, orthogonal_origin);
}

void GradientMesh::remove_edge(Id<HalfEdge> edge)
{
  auto& e = edges[edge];
  auto h = e.handles();
  if (h.has_value() && (!e.twin.has_value() || !edges.contains(e.twin.value())))
  {
    auto [start, end] = h.value();
    handles.remove(start);
    handles.remove(end);
  }
  edges.remove(edge);
}

Id<HalfEdge> GradientMesh::merge_adjacent_siblings(Id<HalfEdge> left,
                                                   Id<HalfEdge> adjacent)
{
  auto orthogonal_edge = eligible_parent(edges[left].next, edges);
  connect(left, edges[adjacent].next);
  auto& f = std::get<Child>(edges[left].kind);
  f.interval.end = edges[adjacent].interval().end;
  remove_edge(adjacent);

  return orthogonal_edge;
}

Id<HalfEdge> GradientMesh::merge_parent_with_child(Id<HalfEdge> parent_or_child)
{
  Id<HalfEdge> parent, child;
  ChildMap child_map;
  if (edges[parent_or_child].leftmost_child.has_value())
  {
    parent = parent_or_child;
    child_map = children(parent, edges);

    if (child_map.size() != 1)
    {
      return parent;
    }
    child = child_map.begin()->second;
  }
  else
  {
    child = parent_or_child;
    parent = eligible_parent(child, edges);
    child_map = children(parent, edges);

    if (child_map.size() != 1)
    {
      return child;
    }
  }

  auto prev = edges[child].prev;
  auto prev_parent = eligible_parent(prev, edges);
  if (prev_parent != prev && edges[prev_parent].next == child)
  {
    connect(prev_parent, parent);
  }
  connect(prev, parent);
  connect(parent, edges[child].next);
  ::connect(parent, edges[child].patch, edges, patches);
  edges[parent].twist = edges[child].twist;
  edges.remove(child);
  edges[parent].leftmost_child = std::nullopt;
  return parent;
}

void GradientMesh::remove_edge_parent_child(Id<HalfEdge> parent)
{
  auto child = edges[parent].leftmost_child;
  if (child.has_value())
  {
    remove_edge_parent_child(child.value());
  }
  remove_handles(parent);
  remove_edge(parent);
}

void GradientMesh::remove_handles(Id<HalfEdge> edge)
{
  auto h = edges[edge].handles();
  if (h.has_value())
  {
    auto [start, end] = h.value();
    handles.remove(start);
    handles.remove(end);
  }
}

void GradientMesh::update_twists(Id<HalfEdge> edge,
                                 std::array<Interpolant, 2> new_twists)
{
  auto& e = edges[edge];
  e.twist = new_twists[0];
  edges[e.next].twist = new_twists[1];
}

void GradientMesh::set_twin(Id<HalfEdge> parent, Id<HalfEdge> twin)
{
  edges[parent].twin = twin;
  auto child_map = children(parent, edges);
  for (auto it = child_map.begin(); it != child_map.end(); ++it)
  {
    edges[it->second].twin = twin;
  }
}

void GradientMesh::share_handles(Id<HalfEdge> source, Id<HalfEdge> twin)
{
  auto h_source = edges[source].handles().value();
  auto h_twin = edges[twin].handles().value();
  if (h_twin[0] != h_source[1] || h_twin[1] != h_source[0])
  {
    handles.remove(h_twin[0]);
    handles.remove(h_twin[1]);
    edges[twin].set_handles({h_source[1], h_source[0]});
  }
}
