#include "gmt/gradient-mesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace gmt
{

using namespace hermite;

// Rotates a 2D vector by 90 degrees clockwise.
static Vector2 rotate(const Vector2& v) { return Vector2(-v.y, v.x); }

GradientMesh::GradientMesh(float side, const std::array<Vector3, 4>& colors)
{
  const auto p = patches.add(Patch{});
  // Tangents colors are initialized to 0 for G^1 continuity.
  // See "Locally Refinable gradient meshes" section 3.2 for reference.
  auto direction = Vector2(1.0f, 0.0f);
  auto start = Vector2(-0.5f, -0.5f) * side;

  // We generate the quad by starting in the top-left corner, advancing in the
  // current direction, and then rotating by 90 degrees,
  // until we loop back to the start.
  auto edge = half_edge(points.add(start),
                        {handles.add({direction}), handles.add({-direction})},
                        colors[0]);

  gmt::connect(edge, p, edges, patches);
  auto top = edge;
  for (int i = 1; i < 4; ++i)
  {
    start += direction * side;
    direction = rotate(direction);
    auto point = points.add(ControlPoint{start});

    auto next =
        half_edge(point, {handles.add({direction}), handles.add({-direction})},
                  colors[i]);

    connect(edge, next);
    edge = next;
  }
  // Close the edge loop
  connect(edge, top);
}

std::pair<Id<Handle>, float> GradientMesh::nearest_handle(
    Vector2 const& position)
{
  float min_distance = std::numeric_limits<float>::infinity();
  Id<Handle> nearest;

  for (auto it = handles.begin(); it != handles.end(); ++it)
  {
    auto [start, end] = endpoints(*it);
    float dist = distance(position, end.coords);

    if (dist < min_distance)
    {
      nearest = handles.get_handle(it);
      min_distance = dist;
    }
  }

  return std::make_pair(nearest, min_distance);
}

std::pair<Id<HalfEdge>, float> GradientMesh::nearest_edge(
    Vector2 const& position)
{
  float min_distance = std::numeric_limits<float>::infinity();
  Id<HalfEdge> nearest;

  // We only consider edges directly adjacent to patches. This way, parent
  // edges that have been split are ignored in favour of their children.
  for (auto const& patch : patches)
  {
    auto top = patch.side;
    auto left = edges[top].prev;
    auto right = edges[top].next;
    auto bottom = edges[right].next;

    float top_dist = distance_to_curve(position, curve_matrix(top));
    float left_dist = distance_to_curve(position, curve_matrix(left));
    float bottom_dist = distance_to_curve(position, curve_matrix(bottom));
    float right_dist = distance_to_curve(position, curve_matrix(right));

    if (top_dist < min_distance)
    {
      nearest = top;
      min_distance = top_dist;
    }

    if (left_dist < min_distance)
    {
      nearest = left;
      min_distance = left_dist;
    }

    if (right_dist < min_distance)
    {
      nearest = right;
      min_distance = right_dist;
    }

    if (bottom_dist < min_distance)
    {
      nearest = bottom;
      min_distance = bottom_dist;
    }
  }

  return std::make_pair(nearest, min_distance);
}

std::pair<Id<ControlPoint>, float> GradientMesh::nearest_point(
    Vector2 const& position)
{
  float min_distance = std::numeric_limits<float>::infinity();
  Id<ControlPoint> nearest;

  for (auto it = points.begin(); it != points.end(); ++it)
  {
    float dist = distance(position, it->coords);

    if (dist < min_distance)
    {
      nearest = points.get_handle(it);
      min_distance = dist;
    }
  }

  return std::make_pair(nearest, min_distance);
}

Vector2 GradientMesh::get_position(Id<ControlPoint> point)
{
  return points[point].coords;
}

void GradientMesh::move_point(Id<ControlPoint> point, Vector2 const& delta)
{
  points[point].coords += delta;
}
Vector2 GradientMesh::get_position(Id<Handle> handle)
{
  auto [start, end] = endpoints(handles[handle]);
  return end.coords;
}

void GradientMesh::move_handle(Id<Handle> handle, Vector2 const& delta)
{
  handles[handle].tangent.coords += 3.0f * delta;
}

void GradientMesh::set_position(Id<HalfEdge> edge, Vector2 const& position)
{
  auto [v1, v2] = get_control_points(edge);
  auto mid_point = 0.5 * (points[v1].coords + points[v2].coords);
  auto v1_offset = points[v1].coords - mid_point;
  auto v2_offset = points[v2].coords - mid_point;

  points[v1].coords = position + v1_offset;
  points[v2].coords = position + v2_offset;
}

Interval GradientMesh::getInterval(Id<HalfEdge> edge) const
{
  return edges[edge].interval();
}

bool GradientMesh::isChild(Id<HalfEdge> edge) const
{
  return !std::get_if<Parent>(&edges[edge].kind);
}

Id<HalfEdge> GradientMesh::next(Id<HalfEdge> edge) const
{
  return edges[edge].next;
}

Id<HalfEdge> GradientMesh::prev(Id<HalfEdge> edge) const
{
  return edges[edge].prev;
}

void GradientMesh::set_color(Id<HalfEdge> edge, Vector3 color)
{
  auto& e = edges[edge];
  e.color = color;
  if (auto h = e.handles())
  {
    handles[h.value()[0]].tangent.color = Vector3();
    handles[h.value()[1]].tangent.color = Vector3();
  }
}

void GradientMesh::set_color_vertex(Id<HalfEdge> edge, Vector3 color)
{
  int cnt = 0;
  auto he = edge;
  while (cnt < 4)
  {
    auto parent_f = std::get_if<Parent>(&edges[he].kind);
    if (parent_f)
    {
      if (edges[he].leftmost_child.has_value())
      {
        set_color(edges[he].leftmost_child.value(), color);
      }
    }
    set_color(he, color);

    if (edges[he].twin.has_value())
    {
      he = edges[edges[he].twin.value()].next;
    }
    else
    {
      break;
    }
    cnt++;
  }
}

int GradientMesh::handle_count() const { return handles.size(); }

int GradientMesh::edge_count() const { return edges.size(); }

int GradientMesh::patch_count() const { return patches.size(); }

int GradientMesh::point_count() const { return points.size(); }

std::vector<PatchRenderData> GradientMesh::render_data() const
{
  std::vector<PatchRenderData> ret;
  ret.reserve(patches.size());

  for (const auto& patch : patches)
    ret.emplace_back(patch_matrix(patch.side).transposed(), boundary(patch));

  return ret;
}

std::vector<PatchMatrix> GradientMesh::patch_data() const
{
  std::vector<PatchMatrix> ret;
  ret.reserve(patches.size());

  for (const auto& patch : patches) ret.emplace_back(patch_matrix(patch.side));

  return ret;
}

std::vector<CurveMatrix> GradientMesh::curve_data() const
{
  std::vector<CurveMatrix> ret;
  ret.reserve(patches.size() * 4);
  for (const auto& patch : patches)
  {
    auto top = patch.side;
    auto left = edges[top].prev;
    auto right = edges[top].next;
    auto bottom = edges[right].next;

    ret.push_back(curve_matrix(top));
    ret.push_back(curve_matrix(left));
    ret.push_back(curve_matrix(bottom));
    ret.push_back(curve_matrix(right));
  }
  return ret;
}

std::vector<Interpolant> GradientMesh::point_data() const
{
  std::vector<Interpolant> ret;
  ret.reserve(points.size());
  for (const auto& point : points)
  {
    auto color = Vector3(0.0f, 0.0f, 0.0f);
    ret.emplace_back(point.coords, color);
  }
  return ret;
}

std::vector<std::array<Interpolant, 2>> GradientMesh::handle_data() const
{
  std::vector<std::array<Interpolant, 2>> ret;
  ret.reserve(handles.size());
  for (const auto& handle : handles)
  {
    ret.push_back(endpoints(handle));
  }
  return ret;
}

Id<HalfEdge> GradientMesh::half_edge(Id<ControlPoint> origin,
                                     std::array<Id<Handle>, 2> edge_handles,
                                     Vector3 color, Interpolant twist,
                                     std::optional<Id<HalfEdge>> twin)
{
  auto edge = edges.add(HalfEdge(origin, edge_handles, color, twist, twin));
  points[origin].edge = edge;
  handles[edge_handles[0]].edge = edge;
  return edge;
}

Id<HalfEdge> GradientMesh::half_edge(
    Id<HalfEdge> parent, Interval interval,
    std::optional<std::array<Id<Handle>, 2>> edge_handles, Vector3 color,
    Interpolant twist, std::optional<Id<HalfEdge>> twin)
{
  auto edge =
      edges.add(HalfEdge(parent, interval, edge_handles, color, twist, twin));
  if (edge_handles)
  {
    handles[edge_handles.value()[0]].edge = edge;
  }
  return edge;
}

void GradientMesh::connect(Id<HalfEdge> a, Id<HalfEdge> b)
{
  edges[a].next = b;
  edges[b].prev = a;
  edges[b].patch = edges[a].patch;
  if (auto h = edges[a].handles())
  {
    auto end = h.value()[1];
    handles[end].edge = b;
  }
}

std::array<Interpolant, 4> GradientMesh::edge_tangents(Id<HalfEdge> edge) const
{
  const auto& e = edges[edge];
  // Origin coordinates
  Vector2 point;
  // Derivatives
  Interpolant start, end;
  visit(
      [&](const Parent& a)
      {
        // Base case: data is directly computable from edge.
        point = points[a.origin].coords;
        start = handles[a.handles[0]].tangent;
        end = -handles[a.handles[1]].tangent;
      },
      [&](const Child& f)
      {
        // Recursive case: data must be computed from parent edge.
        auto mat = curve_matrix(f.parent);
        auto t = f.interval.start;
        point = interpolate(mat, t).coords;
        if (f.handles)
        {
          auto [start_handle, end_handle] = f.handles.value();
          start = handles[start_handle].tangent;
          end = -handles[end_handle].tangent;
        }
        else
        {
          // Half-edge must be parallel to parent
          assert(!is_empty(f.interval));
          auto scale = length(f.interval);
          start = scale * parallel_derivative(mat, t);
          end = scale * parallel_derivative(mat, f.interval.end);
        }
      },
      e);
  auto origin = Interpolant(point, e.color);
  return {origin, e.twist, start, end};
}

PatchMatrix GradientMesh::patch_matrix(Id<HalfEdge> top) const
{
  const auto right = edges[top].next;
  const auto bottom = edges[right].next;
  const auto left = edges[bottom].next;
  assert(top == edges[left].next);

  auto [m0, m0uv, m0v, m1v] = edge_tangents(top);
  auto [m1, m1uv, m1u, m3u] = edge_tangents(right);
  auto [m3, m3uv, m3v, m2v] = edge_tangents(bottom);
  auto [m2, m2uv, m2u, m0u] = edge_tangents(left);

  return {{m0, m0v, m1v, m1,       //
           -m0u, m0uv, -m1uv, m1u, //
           -m2u, -m2uv, m3uv, m3u, //
           m2, -m2v, -m3v, m3}};
}

CurveMatrix GradientMesh::curve_matrix(Id<HalfEdge> edge) const
{
  auto [m0, twist, m0v, m1v] = edge_tangents(edge);
  auto m1 = edge_tangents(edges[edge].next)[0];
  return {{m0, m0v, m1v, m1}};
}

EdgeBoundary GradientMesh::boundary(Id<HalfEdge> edge) const
{
  using Segment = EdgeBoundary::Segment;
  std::vector<Segment> segments;

  auto interval = edges[edge].relative_interval();
  float start = interval.start;
  for (auto [end, unused] : abutting_edges(edge, edges))
  {
    segments.push_back(Segment{{start, end}});
    start = end;
  }
  segments.push_back(Segment{{start, interval.end}});

  auto num_siblings = siblings_abutting_twin(edge, edges);

  auto parent = eligible_parent(edge, edges);
  return EdgeBoundary{curve_matrix(parent), segments, num_siblings};
}

PatchBoundary GradientMesh::boundary(Patch const& patch) const
{
  Id<HalfEdge> top = patch.side;
  Id<HalfEdge> left = edges[top].prev;
  Id<HalfEdge> bottom = edges[left].prev;
  Id<HalfEdge> right = edges[bottom].prev;

  return PatchBoundary{boundary(top), boundary(left), boundary(bottom),
                       boundary(right)};
}

std::array<Interpolant, 2> GradientMesh::endpoints(const Handle& handle) const
{
  auto end = handle.tangent;
  auto start = edge_tangents(handle.edge)[0].coords;
  end.coords = start + end.coords / 3.0f;
  return {start, end};
}

CurveMatrix GradientMesh::split_curve(Id<HalfEdge> edge, float t) const
{
  auto mat = patch_matrix(edge);
  auto p0 = interpolate(mat, 0.0f, t);
  auto m0 = orthogonal_derivative(mat, 0.0f, t);
  auto m1 = orthogonal_derivative(mat, 1.0f, t);
  auto p1 = interpolate(mat, 1.0f, t);
  return {{p0, m0, m1, p1}};
}

std::set<float> GradientMesh::snap_points(Id<HalfEdge> edge) const
{
  std::set<float> ret;
  for (auto [point, e] : abutting_edges(edge, edges, true))
  {
    ret.insert(point);
  }
  for (auto [point, e] : abutting_edges(opposite(edge, edges), edges, true))
  {
    ret.insert(1.0f - point);
  }
  return ret;
}

void GradientMesh::read_curve_matrix(HalfEdge& edge, CurveMatrix const& matrix)
{
  // Set common edge properties.
  edge.color = matrix[0].color;
  edge.twist = matrix[1];

  // Set parent and child specific properties.
  visit(
      [&](Parent const& parent)
      {
        points[parent.origin].coords = matrix[0].coords;
        handles[parent.handles[0]].tangent = matrix[2];
        handles[parent.handles[1]].tangent = -matrix[3];
      },
      [&](Child& child)
      {
        std::array<Id<Handle>, 2> child_handles;

        if (child.handles)
        {
          child_handles = child.handles.value();
        }
        else
        {
          child_handles[0] = handles.add(Handle({}));
          child_handles[1] = handles.add(Handle({}));
          child.handles = child_handles;
        }

        handles[child_handles[0]].tangent = matrix[2];
        handles[child_handles[1]].tangent = -matrix[3];
      },
      edge);
}

void GradientMesh::read_patch_matrix(Patch const& patch,
                                     hermite::PatchMatrix const& matrix)
{
  auto top = patch.side;
  auto right = edges[top].next;
  auto bottom = edges[right].next;
  auto left = edges[bottom].next;
  assert(top == edges[left].next);

  // Following the patch_matrix() function, the target curve matrices are:
  // top:    {m0, m0uv, m0v, m1v}
  // right:  {m1, m1uv, m1u, m3u}
  // bottom: {m3, m3uv, m3v, m2v}
  // left    {m2, m2uv, m2u, m0u}
  //
  // The given patch matrix is:
  // |   m0   m0v   m1v  m1 |
  // | -m0u  m0uv -m1uv m1u |
  // | -m2u -m2uv  m3uv m3u |
  // |   m2  -m2v  -m3v  m3 |

  // Reconstruct curve matrices from the patch matrix.
  CurveMatrix top_curve{
      {matrix(0, 0), matrix(1, 1), matrix(0, 1), matrix(0, 2)}};
  CurveMatrix right_curve{
      {matrix(0, 3), -matrix(1, 2), matrix(1, 3), matrix(2, 3)}};
  CurveMatrix bottom_curve{
      {matrix(3, 3), matrix(2, 2), -matrix(3, 2), -matrix(3, 1)}};
  CurveMatrix left_curve{
      {matrix(3, 0), -matrix(2, 1), -matrix(2, 0), -matrix(1, 0)}};

  read_curve_matrix(edges[top], top_curve);
  read_curve_matrix(edges[right], right_curve);
  read_curve_matrix(edges[bottom], bottom_curve);
  read_curve_matrix(edges[left], left_curve);
}

void GradientMesh::read_patch_data(std::vector<PatchMatrix> const& data)
{
  assert(patches.size() == data.size());

  size_t index = 0;
  for (auto const& patch : patches)
  {
    read_patch_matrix(patch, data[index]);
    ++index;
  }
}

} // namespace gmt
