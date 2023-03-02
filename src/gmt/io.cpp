#include <fstream>
#include <sstream>

#include "gmt/gradient-mesh.hpp"

namespace gmt
{

using namespace hermite;

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

    tokens >> points[point].coords.x;
    tokens >> points[point].coords.y;

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

    tokens >> handles[handle].tangent.coords.x;
    tokens >> handles[handle].tangent.coords.y;

    tokens >> handles[handle].tangent.color.r;
    tokens >> handles[handle].tangent.color.g;
    tokens >> handles[handle].tangent.color.b;
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
    tokens >> twist.coords.x >> twist.coords.y;
    tokens >> twist.color.r >> twist.color.g >> twist.color.b;

    Vector3 color;
    tokens >> color.r >> color.g >> color.b;

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

      edge =
          half_edge(parent, child_interval, child_handles, color, twist, twin);
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
    output << point.edge.id << '\n';
  }
}

void GradientMesh::write_handles(std::ostream &output) const
{
  for (auto const &handle : handles)
  {
    output << handle.edge.id << ' ';
    output << handle.tangent.coords.x << ' ';
    output << handle.tangent.coords.y << ' ';
    output << handle.tangent.color.r << ' ';
    output << handle.tangent.color.g << ' ';
    output << handle.tangent.color.b << '\n';
  }
}

void GradientMesh::write_patches(std::ostream &output) const
{
  for (auto const &patch : patches) output << patch.side.id << '\n';
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
      output << edge.handles().value()[0].id << ' '; // 10
      output << edge.handles().value()[1].id << ' '; // 11
    }
    else
    {
      output << -1 << ' '; // 10
      output << -1 << ' '; // 11
    }

    if (edge.twin.has_value())
      output << edge.twin.value().id << ' '; // 12
    else
      output << -1 << ' '; // 12

    output << edge.prev.id << ' ';  // 13
    output << edge.next.id << ' ';  // 14
    output << edge.patch.id << ' '; // 15

    if (edge.leftmost_child.has_value())
      output << edge.leftmost_child.value().id << ' '; // 16
    else
      output << -1 << ' '; // 16

    auto child_f = std::get_if<Child>(&edge.kind);
    auto parent_f = std::get_if<Parent>(&edge.kind);
    if (child_f)
    {
      output << 1 << ' ';                  // 17
      output << child_f->parent.id << ' '; // 18

      if (child_f->handles.has_value())
      {
        output << child_f->handles.value()[0].id << ' '; // 19
        output << child_f->handles.value()[1].id << ' '; // 20
      }
      else
      {
        output << -1 << ' '; // 19
        output << -1 << ' '; // 20
      }

      output << child_f->interval.start << ' '; // 21
      output << child_f->interval.end << ' ';   // 22
    }
    else
    {
      output << 0 << ' ';                       // 17
      output << parent_f->handles[0].id << ' '; // 18
      output << parent_f->handles[1].id << ' '; // 19
      output << parent_f->origin.id << ' ';     // 20
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
