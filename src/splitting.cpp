#include <algorithm>
#include <cassert>
#include <cmath>

#include "gradient-mesh.hpp"
#include "interval.hpp"
#include "snapping.hpp"

using namespace hermite;

std::optional<Id<HalfEdge>> GradientMesh::is_parent_t_junction(
    Id<HalfEdge> edge)
{
  auto child_f = std::get_if<Child>(&edges[edge].kind);
  if (child_f && is_full(child_f->interval))
  {
    auto parent_f = std::get_if<Child>(&edges[child_f->parent].kind);
    if (parent_f)
    {
      return child_f->parent;
    }
  }

  return std::nullopt;
}

std::vector<Id<HalfEdge>> GradientMesh::find_t_junctions(Id<HalfEdge> edge)
{
  std::vector<Id<HalfEdge>> t_edges;

  auto parent_f = is_parent_t_junction(edge);
  if (parent_f.has_value())
  {
    t_edges.push_back(parent_f.value());
  }
  else
  {
    auto children = siblings(edge, edges);

    for (auto it = children.begin(); it != children.end(); ++it)
    {
      Id<HalfEdge> current = it->second;
      auto interval = edges[current].interval();
      if (is_empty(interval))
      {
        t_edges.push_back(current);
      }
    }
  }

  return t_edges;
}

void GradientMesh::remove_t_junctions()
{
  for (auto it = patches.begin(); it != patches.end(); ++it)
  {
    auto edge = it->side;
    for (size_t i = 0; i < 4; ++i)
    {
      std::vector<Id<HalfEdge>> t_edges = find_t_junctions(edge);
      for (auto current : t_edges)
      {
        auto child_f = std::get_if<Child>(&edges[current].kind);
        if (child_f)
        {
          float t = edges[current].interval().start;
          propagate_split(child_f->parent, 1.0 - t);

          // just reset iterator what are you gonna do stab me?
          it = patches.begin();
          break;
        }
      }

      edge = edges[edge].next;
    }
  }
}

Id<HalfEdge> GradientMesh::split_global(Id<HalfEdge> selected_edge, float t,
                                        float t2)
{
  auto top = selected_edge;
  auto right = edges[selected_edge].next;
  auto bottom = edges[right].next;
  auto left = edges[bottom].next;

  propagate_split(top, 1.0 - t);
  propagate_split(bottom, t);
  propagate_split(right, 1.0 - t2);
  propagate_split(left, t2);

  return subdivide(selected_edge, t, t2);
}

void GradientMesh::propagate_split(Id<HalfEdge> top, float t)
{
  auto twin = edges[top].twin;
  auto prev_edge = top;
  float next_t;

  t = calculate_next_t(top, t);

  while (twin.has_value())
  {
    auto edge_to_split = twin.value();
    auto opp = opposite(edge_to_split, edges);
    auto next_twin = edges[opp].twin;

    if (edges[edge_to_split].leftmost_child.has_value())
    {
      auto child = find_child_contains_t(edge_to_split, t);
      t = find_t_on_child(edge_to_split, t);
      opp = opposite(child, edges);
      next_twin = edges[opp].twin;

      next_t = calculate_next_t(opp, t);

      split(child, t);
    }
    else
    {
      next_t = calculate_next_t(opp, t);

      split(edge_to_split, t);
    }

    twin = next_twin;
    opp = prev_edge;
    t = next_t;
  }
}

float GradientMesh::calculate_t_on_parent(Id<HalfEdge> child, float t)
{
  auto interval = edges[child].interval();
  float t0 = interval.start;
  float t1 = interval.end;
  t = t * t0 + (1.0 - t) * t1;
  return 1.0 - t;
}

float GradientMesh::calculate_next_t(Id<HalfEdge> opp, float t)
{
  float next_t = t;
  auto opp_f = std::get_if<Child>(&edges[opp].kind);
  if (opp_f && !is_empty(opp_f->interval))
  {
    next_t = calculate_t_on_parent(opp, t);
  }
  return next_t;
}

Id<HalfEdge> GradientMesh::find_child_contains_t(Id<HalfEdge> parent, float t)
{
  auto siblings = children(parent, edges);

  Id<HalfEdge> child_with_interval = parent;
  for (auto it = siblings.begin(); it != siblings.end(); ++it)
  {
    auto child = it->second;
    Interval interval = edges[child].relative_interval();

    if (contains(interval, t))
    {
      child_with_interval = child;
      return child;
    }
  }
  return child_with_interval;
}

float GradientMesh::find_t_on_child(Id<HalfEdge> parent, float t)
{
  auto child = find_child_contains_t(parent, t);
  Interval interval = edges[child].relative_interval();
  return relative_child_position(interval, t);
}

Id<HalfEdge> GradientMesh::subdivide(Id<HalfEdge> selected_edge, float t,
                                     float t2)
{
  auto next = edges[selected_edge].next;
  auto prev = edges[selected_edge].prev;
  selected_edge = split(selected_edge, t);

  next = split(next, t2);
  prev = split(prev, 1.0 - t2);

  return selected_edge;
}

Id<HalfEdge> GradientMesh::split(Id<HalfEdge> top, float t)
{
  auto right = edges[top].next;
  auto bottom = edges[right].next;
  auto left = edges[bottom].next;
  auto matrix = patch_matrix(top);

  auto top_handle = handles.add(orthogonal_derivative(matrix, 0.0f, t));
  auto bottom_handle = handles.add(-orthogonal_derivative(matrix, 1.0f, t));

  auto [top_left, middle_left, top_right] = junction(
      top, t, mixed_derivative(matrix, 0.0f, t), {top_handle, bottom_handle});

  auto [bottom_right, middle_right, bottom_left] =
      junction(bottom, 1.0f - t, mixed_derivative(matrix, 1.0f, t),
               {bottom_handle, top_handle});

  // Link references in the left patch
  auto left_patch = edges[top].patch;
  ::connect(top_left, left_patch, edges, patches);
  connect(top_left, middle_left);
  connect(middle_left, bottom_left);
  connect(bottom_left, left);
  connect(left, top_left);

  // Link references in the right patch
  auto right_patch = patches.add(Patch{});
  ::connect(top_right, right_patch, edges, patches);
  connect(top_right, right);
  connect(right, bottom_right);
  connect(bottom_right, middle_right);
  connect(middle_right, top_right);

  scale_twists(top_left, t);
  scale_twists(bottom_left, t);
  scale_twists(top_right, 1.0f - t);
  scale_twists(bottom_right, 1.0f - t);

  connect_twins(middle_left, middle_right, edges);

  snap_to_twin(top_right);
  snap_to_twin(bottom_left);

  return top_left;
}

std::array<Id<HalfEdge>, 3> GradientMesh::junction(
    Id<HalfEdge> edge, float t, Interpolant twist,
    std::array<Id<Handle>, 2> ortho_handles)
{
  if (edges[edge].twin.has_value())
  {
    return t_junction(edge, t, twist, ortho_handles);
  }
  else
  {
    return connected_junction(edge, t, twist, ortho_handles);
  }
}

std::array<Id<HalfEdge>, 3> GradientMesh::connected_junction(
    Id<HalfEdge> left, float t, Interpolant twist,
    std::array<Id<Handle>, 2> ortho_handles)
{
  // A half-edge with no twin must be attached.
  auto& e = edges[left];
  auto& a = std::get<Parent>(e.kind);
  const auto [start_handle, end_handle] = a.handles;
  auto twin = e.twin;

  auto matrix = curve_matrix(left);
  auto mid = interpolate(matrix, t);
  auto deriv = parallel_derivative(matrix, t);

  a.handles[1] = handles.add(-deriv);
  auto mid_point = points.add(mid.coords);
  auto right = half_edge(mid_point, {handles.add(deriv), end_handle}, mid.color,
                         twist, twin);
  auto middle = half_edge(mid_point, ortho_handles, mid.color, -twist);

  scale_tangents(left, t);
  scale_tangents(right, 1.0f - t);

  return {left, middle, right};
}

void GradientMesh::scale_tangents(Id<HalfEdge> edge, float t)
{
  auto& e = edges[edge];
  if (auto h = e.handles())
  {
    auto [start, end] = h.value();
    handles[start].tangent *= t;
    handles[end].tangent *= t;
  }
}

void GradientMesh::scale_twists(Id<HalfEdge> edge, float t)
{
  edges[edge].twist *= t;
  edges[edges[edge].next].twist *= t;
}

std::array<Id<HalfEdge>, 3> GradientMesh::t_junction(
    Id<HalfEdge> edge, float t, Interpolant twist,
    std::array<Id<Handle>, 2> ortho_handles)
{
  auto color = interpolate(curve_matrix(edge), t).color;

  const auto e = edges[edge];
  const auto parent = eligible_parent(edge, edges);
  const auto interval = e.relative_interval();
  auto mid = interval(t);

  auto left = child_from(edge, mid);
  auto right = half_edge(parent, Interval(mid, interval.end), std::nullopt,
                         color, twist, e.twin);
  auto middle =
      half_edge(parent, Interval(mid, mid), ortho_handles, color, -twist);

  return {left, middle, right};
}

Id<HalfEdge> GradientMesh::child_from(Id<HalfEdge> edge, float end_point)
{
  auto& e = edges[edge];
  if (e.is_parentable())
  {
    auto ret = half_edge(edge, Interval(0.0f, end_point), std::nullopt, e.color,
                         e.twist, e.twin);
    edges[edge].leftmost_child = ret;
    return ret;
  }
  else
  {
    auto& f = std::get<Child>(e.kind);
    f.interval.end = end_point;
    return edge;
  }
}

void GradientMesh::to_parent(Id<HalfEdge> edge, Id<ControlPoint> origin)
{
  auto& e = edges[edge];
  // Edge must have handles to be attached.
  e.kind = Parent{origin, e.handles().value()};
}

void GradientMesh::snap_to_twin(Id<HalfEdge> edge)
{
  const auto& e = edges[edge];
  if (!e.twin.has_value() || e.is_parentable())
  {
    return;
  }
  auto parent = eligible_parent(edge, edges);
  auto siblings = children(parent, edges);
  auto twin_parent = e.twin.value();
  auto twins = children(twin_parent, edges);

  auto interval = e.relative_interval();

  auto edge_it = find(edge, siblings);
  auto twin_it = find_approx(twins, converse(interval).end);
  assert(edge_it != siblings.end() && edge_it->second == edge);
  if (twin_it == twins.end())
  {
    return;
  }
  auto twin = twin_it->second;

  // Don't snap when there is an offset between the edge and its twin.
  auto offset = (edge_it->first - (1.0f - twin_it->first));
  if (std::abs(offset) > MAXIMUM_TWIN_SNAPPING_OFFSET)
  {
    return;
  }

  auto t = interval.start;
  auto matrix = curve_matrix(parent);
  auto mid = interpolate(matrix, t);
  auto parallel = parallel_derivative(matrix, t);
  auto twist = edges[edge].twist;

  auto mid_point = points.add(mid.coords);
  auto left_handle = handles.add(-parallel);
  auto right_handle = handles.add(parallel);

  auto bottom_edge = edges[edges[edge].prev].twin.value();
  auto top_edge = edges[edges[twin].prev].twin.value();

  // Split the parent edge
  auto [left, right] = parent_junction(parent, t, mid_point, mid.color, twist,
                                       {left_handle, right_handle});
  edges[left].next = bottom_edge;
  to_parent(bottom_edge, mid_point);

  auto [twin_right, twin_left] =
      parent_junction(twin_parent, 1.0f - t, mid_point, mid.color, twist,
                      {right_handle, left_handle});
  edges[twin_right].next = top_edge;
  to_parent(top_edge, mid_point);

  // Revert the effect of scaling handles twice
  // due to the two calls to connected_junction
  scale_tangents(left, 1.0f / t);
  scale_tangents(right, 1.0f / (1.0f - t));

  connect_twins(left, twin_left, edges);
  connect_twins(right, twin_right, edges);

  reparametrize(siblings.begin(), edge_it, left, Interval(0.0f, t));
  reparametrize(edge_it, siblings.end(), right, Interval(t, 1.0f));

  reparametrize(twins.begin(), twin_it, twin_right, Interval(0.0f, 1.0f - t));
  reparametrize(twin_it, twins.end(), twin_left, Interval(1.0f - t, 1.0f));
}

std::array<Id<HalfEdge>, 2> GradientMesh::parent_junction(
    Id<HalfEdge> left, float t, Id<ControlPoint> mid_point, Vector3 color,
    Interpolant twist, std::array<Id<Handle>, 2> parallel_handles)
{
  auto& e = edges[left];
  const auto [start_handle, end_handle] = e.handles().value();
  const auto twin = e.twin;
  const auto next = e.next;

  e.twist *= t;
  e.set_handles({start_handle, parallel_handles[0]});
  auto right = half_edge(mid_point, {parallel_handles[1], end_handle}, color,
                         twist * (1.0f - t), twin);

  // Only needed to find edge endpoint.
  edges[right].next = next;

  scale_tangents(left, t);
  scale_tangents(right, 1.0f - t);

  return {left, right};
}

void GradientMesh::reparametrize(ChildIterator begin, ChildIterator end,
                                 Id<HalfEdge> parent, Interval new_interval)
{
  edges[parent].leftmost_child = begin->second;
  auto twin = edges[parent].twin;
  for (auto it = begin; it != end; ++it)
  {
    auto& e = edges[it->second];
    e.twin = twin;

    auto& f = std::get<Child>(e.kind);
    auto f_original_parent = f.parent;
    f.interval = ::reparametrize(f.interval, new_interval);
    f.parent = parent;

    auto next = eligible_parent(e.next, edges);
    auto next_f = std::get_if<Child>(&edges[next].kind);
    // Change orthogonal edges only when they share a parent with the current
    // child.
    if (next_f && is_empty(next_f->interval) &&
        next_f->parent == f_original_parent)
    {
      next_f->interval = ::reparametrize(next_f->interval, new_interval);
      next_f->parent = parent;
    }
  }
}
