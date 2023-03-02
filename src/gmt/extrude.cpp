#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "gmt/gradient-mesh.hpp"

namespace gmt
{

using namespace hermite;

Id<ControlPoint> GradientMesh::get_origin(Id<HalfEdge> edge)
{
  Id<ControlPoint> origin;

  auto parent_f = std::get_if<Parent>(&edges[edge].kind);
  if (parent_f)
  {
    origin = parent_f->origin;
  }
  else
  {
    auto child_f = std::get_if<Child>(&edges[edge].kind);
    if (child_f)
    {
      origin = get_origin(child_f->parent);
    }
  }
  return origin;
}

std::array<Id<ControlPoint>, 2> GradientMesh::get_control_points(
    Id<HalfEdge> edge)
{
  auto next = edges[edge].next;

  return {get_origin(edge), get_origin(next)};
}

Id<Handle> GradientMesh::get_handle_by_number(Id<HalfEdge> edge, int index)
{
  Id<Handle> handle;
  if (edges[edge].handles().has_value())
  {
    handle = edges[edge].handles().value()[index];
  }
  else
  {
    auto child_f = std::get_if<Child>(&edges[edge].kind);
    if (child_f)
    {
      handle = edges[child_f->parent].handles().value()[index];
    }
  }
  return handle;
}

Id<HalfEdge> GradientMesh::extrude(Id<HalfEdge> edge)
{
  if (edges[edge].twin.has_value())
  {
    return edge;
  }

  auto next = edges[edge].next;
  auto prev = edges[edge].prev;
  auto [v1, v2] = get_control_points(edge);

  auto new_patch = patches.add(Patch{});
  auto v1_prime = points.add(points[v1].coords);
  auto v2_prime = points.add(points[v2].coords);

  auto twin_handles = edges[edge].handles().value();
  std::reverse(twin_handles.begin(), twin_handles.end());

  auto next_handle = get_handle_by_number(next, 0);
  auto prev_handle = get_handle_by_number(prev, 1);

  auto next_twin_handle = handles.add(Handle(-handles[prev_handle].tangent));
  auto next_twin_handle_2 = handles.add(Handle(handles[prev_handle].tangent));

  auto prev_twin_handle = handles.add(Handle(handles[next_handle].tangent));
  auto prev_twin_handle_2 = handles.add(Handle(-handles[next_handle].tangent));

  auto opp_handle = handles.add(Handle(handles[twin_handles[1]]));
  auto opp_handle_2 = handles.add(Handle(handles[twin_handles[0]]));

  Vector3 color_v2 = edges[next].color;
  Vector3 color_v1 = edges[edge].color;
  auto twist = edges[edge].twist;

  auto twin = half_edge(v2, twin_handles, color_v2, twist, edge);
  auto twin_next =
      half_edge(v1, {next_twin_handle, next_twin_handle_2}, color_v1, twist);
  auto twin_prev = half_edge(v2_prime, {prev_twin_handle, prev_twin_handle_2},
                             color_v2, twist);
  auto twin_opp =
      half_edge(v1_prime, {opp_handle, opp_handle_2}, color_v1, twist);

  gmt::connect(twin, new_patch, edges, patches);
  gmt::connect(twin_next, new_patch, edges, patches);
  gmt::connect(twin_prev, new_patch, edges, patches);
  gmt::connect(twin_opp, new_patch, edges, patches);

  connect_twins(edge, twin, edges);

  connect(twin, twin_next);
  connect(twin_next, twin_opp);
  connect(twin_opp, twin_prev);
  connect(twin_prev, twin);

  return twin_opp;
}

} // namespace gmt
