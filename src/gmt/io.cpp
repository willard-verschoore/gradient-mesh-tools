#include <fstream>
#include <sstream>

#include "gmt/gradient-mesh.hpp"

namespace gmt
{

using namespace hermite;

size_t GradientMesh::get_index(Id<Handle> handle) const
{
  return handles.find(handle) - handles.begin();
}

size_t GradientMesh::get_index(Id<HalfEdge> edge) const
{
  return edges.find(edge) - edges.begin();
}

size_t GradientMesh::get_index(Id<Patch> patch) const
{
  return patches.find(patch) - patches.begin();
}

size_t GradientMesh::get_index(Id<ControlPoint> point) const
{
  return points.find(point) - points.begin();
}

// Utility function for reading 2D positions.
static void read_position(std::istream &input, Vector2 &position)
{
  input >> position.x;
  input >> position.y;
}

// Utility function for reading RGB colors. Can Clamp the channels to [0, 1].
static void read_color(std::istream &input, Vector3 &color, bool clamp)
{
  input >> color.r;
  input >> color.g;
  input >> color.b;

  if (clamp)
  {
    color.r = std::clamp(color.r, 0.0f, 1.0f);
    color.g = std::clamp(color.g, 0.0f, 1.0f);
    color.b = std::clamp(color.b, 0.0f, 1.0f);
  }
}

static void read_header(std::istream &input, int &num_points, int &num_handles,
                        int &num_patches, int &num_edges)
{
  std::string line;
  getline(input, line); // Line 1: HEMESH.
  getline(input, line); // Line 2: points handles patches edges.
  std::istringstream tokens(line);

  tokens >> num_points >> num_handles >> num_patches >> num_edges;
}

void GradientMesh::read_points(std::istream &input, int num_points)
{
  std::string line;

  for (int i = 0; i < num_points; ++i)
  {
    getline(input, line);
    if (line.empty()) continue;
    std::istringstream tokens(line);

    Id<ControlPoint> point = points.add(ControlPoint({}));

    read_position(tokens, points[point].coords);

    Id<HalfEdge> edge{0, 0};
    tokens >> edge.id;
    points[point].edge = edge;
  }
}

void GradientMesh::read_handles(std::istream &input, int num_handles)
{
  std::string line;

  for (int i = 0; i < num_handles; ++i)
  {
    getline(input, line);
    if (line.empty()) continue;
    std::istringstream tokens(line);

    Id<Handle> handle = handles.add(Handle({}));

    Id<HalfEdge> edge{0, 0};
    tokens >> edge.id;
    handles[handle].edge = edge;

    read_position(tokens, handles[handle].tangent.coords);
    read_color(tokens, handles[handle].tangent.color, false);
  }
}

void GradientMesh::read_patches(std::istream &input, int num_patches)
{
  std::string line;

  for (int i = 0; i < num_patches; ++i)
  {
    getline(input, line);
    if (line.empty()) continue;
    std::istringstream tokens(line);

    Id<Patch> patch = patches.add(Patch({}));

    Id<HalfEdge> side{0, 0};
    tokens >> side.id;
    patches[patch].side = side;
  }
}

void GradientMesh::read_edges(std::istream &input, int num_edges)
{
  std::string line;

  for (int i = 0; i < num_edges; ++i)
  {
    getline(input, line);
    if (line.empty()) continue;
    std::istringstream tokens(line);

    Id<HalfEdge> edge;

    // TODO: This is unused. Why? Can it be removed?
    Interval interval{0, 0};
    tokens >> interval.start >> interval.end;

    Interpolant twist;
    read_position(tokens, twist.coords);
    read_color(tokens, twist.color, false);

    Vector3 color;
    read_color(tokens, color, true);

    int handle_1_id, handle_2_id;
    tokens >> handle_1_id >> handle_2_id;
    std::optional<std::array<Id<Handle>, 2>> edge_handles = std::nullopt;
    if (handle_1_id >= 0)
    {
      Id<Handle> handle_1{(uint32_t)handle_1_id, 0};
      Id<Handle> handle_2{(uint32_t)handle_2_id, 0};
      edge_handles = {handle_1, handle_2};
    }

    int twin_id;
    tokens >> twin_id;
    std::optional<Id<HalfEdge>> twin = std::nullopt;
    if (twin_id >= 0) twin = {(uint32_t)twin_id, 0};

    Id<HalfEdge> prev{0, 0}, next{0, 0};
    tokens >> prev.id >> next.id;

    Id<Patch> patch{0, 0};
    tokens >> patch.id;

    int leftmost_child_id;
    tokens >> leftmost_child_id;
    std::optional<Id<HalfEdge>> leftmost_child = std::nullopt;
    if (leftmost_child_id >= 0)
      leftmost_child = {(uint32_t)leftmost_child_id, 0};

    int is_child;
    tokens >> is_child;

    if (is_child)
    {
      Id<HalfEdge> parent{0, 0};
      tokens >> parent.id;

      int handle_1_id, handle_2_id;
      tokens >> handle_1_id >> handle_2_id;
      std::optional<std::array<Id<Handle>, 2>> child_handles = std::nullopt;
      if (handle_1_id >= 0)
      {
        Id<Handle> handle_1{(uint32_t)handle_1_id, 0};
        Id<Handle> handle_2{(uint32_t)handle_2_id, 0};
        child_handles = {handle_1, handle_2};
      }

      Interval child_interval{0, 0};
      tokens >> child_interval.start >> child_interval.end;

      bool recolored = false;
      tokens >> recolored;

      edge = half_edge(parent, child_interval, child_handles, recolored, color,
                       twist, twin);
    }
    else
    {
      Id<Handle> handle_1{0, 0}, handle_2{0, 0};
      tokens >> handle_1.id >> handle_2.id;

      Id<ControlPoint> origin{0, 0};
      tokens >> origin.id;

      edge = half_edge(origin, {handle_1, handle_2}, color, twist, twin);
    }

    edges[edge].next = next;
    edges[edge].prev = prev;
    edges[edge].patch = patch;
    edges[edge].leftmost_child = leftmost_child;
  }
}

void GradientMesh::read_from_file(const std::string &file_name)
{
  std::ifstream file(file_name);
  if (!file.is_open()) return;

  points.clear();
  handles.clear();
  patches.clear();
  edges.clear();

  int num_points, num_handles, num_patches, num_edges;
  read_header(file, num_points, num_handles, num_patches, num_edges);

  read_points(file, num_points);
  read_handles(file, num_handles);
  read_patches(file, num_patches);
  read_edges(file, num_edges);

  file.close();
}

static void write_header(std::ostream &output, int num_points, int num_handles,
                         int num_patches, int num_edges)
{
  output << "HEMESH\n";
  output << num_points << ' ';
  output << num_handles << ' ';
  output << num_patches << ' ';
  output << num_edges << '\n';
}

void GradientMesh::write_points(std::ostream &output) const
{
  for (auto const &point : points)
  {
    output << point.coords.x << ' ';
    output << point.coords.y << ' ';
    output << get_index(point.edge) << '\n';
  }
}

void GradientMesh::write_handles(std::ostream &output) const
{
  for (auto const &handle : handles)
  {
    output << get_index(handle.edge) << ' ';
    output << handle.tangent.coords.x << ' ';
    output << handle.tangent.coords.y << ' ';
    output << handle.tangent.color.r << ' ';
    output << handle.tangent.color.g << ' ';
    output << handle.tangent.color.b << '\n';
  }
}

void GradientMesh::write_patches(std::ostream &output) const
{
  for (auto const &patch : patches) output << get_index(patch.side) << '\n';
}

void GradientMesh::write_edges(std::ostream &output) const
{
  for (auto const &edge : edges)
  {
    output << edge.interval().start << ' '; // 0
    output << edge.interval().end << ' ';   // 1
    output << edge.twist.coords.x << ' ';   // 2
    output << edge.twist.coords.y << ' ';   // 3
    output << edge.twist.color.x << ' ';    // 4
    output << edge.twist.color.y << ' ';    // 5
    output << edge.twist.color.z << ' ';    // 6
    output << edge.color.x << ' ';          // 7
    output << edge.color.y << ' ';          // 8
    output << edge.color.z << ' ';          // 9

    if (edge.handles().has_value())
    {
      output << get_index(edge.handles().value()[0]) << ' '; // 10
      output << get_index(edge.handles().value()[1]) << ' '; // 11
    }
    else
    {
      output << -1 << ' '; // 10
      output << -1 << ' '; // 11
    }

    if (edge.twin.has_value())
      output << get_index(edge.twin.value()) << ' '; // 12
    else
      output << -1 << ' '; // 12

    output << get_index(edge.prev) << ' ';  // 13
    output << get_index(edge.next) << ' ';  // 14
    output << get_index(edge.patch) << ' '; // 15

    if (edge.leftmost_child.has_value())
      output << get_index(edge.leftmost_child.value()) << ' '; // 16
    else
      output << -1 << ' '; // 16

    auto child_f = std::get_if<Child>(&edge.kind);
    auto parent_f = std::get_if<Parent>(&edge.kind);
    if (child_f)
    {
      output << 1 << ' ';                          // 17
      output << get_index(child_f->parent) << ' '; // 18

      if (child_f->handles.has_value())
      {
        output << get_index(child_f->handles.value()[0]) << ' '; // 19
        output << get_index(child_f->handles.value()[1]) << ' '; // 20
      }
      else
      {
        output << -1 << ' '; // 19
        output << -1 << ' '; // 20
      }

      output << child_f->interval.start << ' '; // 21
      output << child_f->interval.end << ' ';   // 22

      output << child_f->recolored << ' '; // 23
    }
    else
    {
      output << 0 << ' ';                               // 17
      output << get_index(parent_f->handles[0]) << ' '; // 18
      output << get_index(parent_f->handles[1]) << ' '; // 19
      output << get_index(parent_f->origin) << ' ';     // 20
    }

    output << '\n';
  }
}

void GradientMesh::write_to_file(std::string const &file_name) const
{
  std::ofstream file(file_name);
  if (!file.is_open()) return;

  write_header(file, points.size(), handles.size(), patches.size(),
               edges.size());

  write_points(file);
  write_handles(file);
  write_patches(file);
  write_edges(file);

  file.close();
}

} // namespace gmt
