#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

#include "gradient-mesh.hpp"

static std::vector<std::string> get_tokens(std::string const &input,
                                           char delimiter)
{
  std::vector<std::string> result;
  std::stringstream stream(input);
  std::string token;

  while (getline(stream, token, delimiter))
    if (!token.empty()) result.push_back(token);

  return result;
}

void read_cgm_handle(std::vector<std::string> const &tokens,
                     SlotMap<Handle> &handles)
{
  Id<HalfEdge> edge = {(uint32_t)std::stoul(tokens[3]), 0};
  Vector2 coord = Vector2(std::stof(tokens[1]), std::stof(tokens[2]));

  auto handle = handles.add(Handle({}));
  handles[handle].edge = edge;
  handles[handle].tangent.coords = coord;

  /*
  // ControlPoint,                                             vCoords,
  vCurveLink,          vValency,
  //gradientMesh->controlPoints.append(ControlPoint(Vector2(values[1].toFloat(),
  values[2].toFloat()), values[3].toInt(), values[4].toInt(),
  // vIsTpoint,                          vTpointConnections,     vLevel, vIndex
  //  (values[5] == "1" ? true : false), values[6].toInt(), values[7].toInt(),
  values[8].toInt()));
  */
}

void read_cgm_point(std::vector<std::string> const &tokens,
                    std::vector<Id<HalfEdge>> &t_stems,
                    SlotMap<ControlPoint> &points)
{
  Id<HalfEdge> edge = {(uint32_t)std::stoul(tokens[3]), 0};
  Vector2 coord = Vector2(std::stof(tokens[1]), std::stof(tokens[2]));

  auto point = points.add(ControlPoint({}));
  points[point].coords = coord;
  points[point].edge = edge;

  bool is_t_point = (tokens[5] == "1");
  if (is_t_point)
  {
    Id<HalfEdge> stem_edge = edge;
    t_stems.push_back(stem_edge);
  }

  /*
  // ControlPoint,                                             vCoords,
  vCurveLink,          vValency,
  //gradientMesh->controlPoints.append(ControlPoint(Vector2(values[1].toFloat(),
  values[2].toFloat()), values[3].toInt(), values[4].toInt(),
  // vIsTpoint,                          vTpointConnections,     vLevel, vIndex
  //  (values[5] == "1" ? true : false), values[6].toInt(), values[7].toInt(),
  values[8].toInt()));
  */
}

bool is_stem(Id<HalfEdge> edge, std::vector<Id<HalfEdge>> const &t_stems)
{
  return std::find(t_stems.begin(), t_stems.end(), edge) != t_stems.end();
}

void GradientMesh::open_from_cgm(std::string const &file_name)
{
  std::ifstream file(file_name);
  if (!file.is_open()) return;

  std::string line;
  std::vector<std::string> tokens;

  points.clear();
  handles.clear();
  edges.clear();
  patches.clear();

  std::set<Id<HalfEdge>> parent_twins;
  std::vector<Id<HalfEdge>> parentless;
  std::vector<Id<HalfEdge>> t_stems;

  uint edge_cnt = 0;

  while (getline(file, line))
  {
    tokens = get_tokens(line, ' ');

    if (tokens[0] == "S")
    {
      // Size of ControlPoints, HalfCurves and Patches
      // num_points = (uint32_t)std::stoul(tokens[1]);
      // num_edges = (uint32_t)std::stoul(tokens[2]);
      // num_patches = (uint32_t)std::stoul(tokens[3]);
    }
    else if (tokens[0] == "H")
    {
      read_cgm_handle(tokens, handles);
    }
    else if (tokens[0] == "C")
    {
      read_cgm_point(tokens, t_stems, points);
    }
    else if (tokens[0] == "E")
    {
      Id<Patch> patch = {(uint32_t)std::stoul(tokens[22]), 0};

      // check if has twin we dont use boundary edges here
      if (patch.id < 65535)
      { // for some reason no patch is MAX_SHORT in pieters code
        Id<ControlPoint> corner = {(uint32_t)std::stoul(tokens[1]), 0};
        Id<HalfEdge> next = {(uint32_t)std::stoul(tokens[19]), 0};
        Id<HalfEdge> prev = {(uint32_t)std::stoul(tokens[20]), 0};

        std::optional<Id<HalfEdge>> twin = std::nullopt;
        if (std::stoi(tokens[21]) >= 0)
        {
          twin = {(uint32_t)std::stoul(tokens[21]), 0};
        }

        Interpolant twist;
        twist.coords = Vector2(std::stof(tokens[4]), std::stof(tokens[5]));
        twist.color = Vector3(std::stof(tokens[15]), std::stof(tokens[16]),
                              std::stof(tokens[17]));

        Vector3 color = Vector3(std::stof(tokens[6]), std::stof(tokens[7]),
                                std::stof(tokens[8]));

        Id<Handle> handleTail = {(uint32_t)std::stoul(tokens[2]), 0};
        Id<Handle> handleHead = {(uint32_t)std::stoul(tokens[3]), 0};
        std::array<Id<Handle>, 2> edge_handles = {handleTail, handleHead};

        Id<HalfEdge> edge;
        int level = std::stoi(tokens[23]);
        if (level > 0)
        {
          Interval interval = {0.0f, 1.0 / pow(2.0, double(level))};
          // Interval interval = {1.0 / pow(2.0, double(level)), 0.0};

          if (twin.has_value())
          {
            if (is_stem(twin.value(), t_stems))
            {
              edge = half_edge(corner, edge_handles, color, twist, twin);
            }
            else
            {
              if (!is_stem({edge_cnt, 0}, t_stems))
              {
                edge = half_edge(Id<HalfEdge>(), interval, edge_handles, color,
                                 twist, twin);

                parentless.push_back(edge);
                parent_twins.insert(twin.value());
              }
              else
              {
                edge = half_edge(Id<HalfEdge>(), {0.5, 0.5}, edge_handles,
                                 color, twist, twin);
              }
            }
          }
          else
          {
            // boundary edge at higher level
            edge = half_edge(corner, edge_handles, color, twist, twin);
          }
        }
        else
        {
          edge = half_edge(corner, edge_handles, color, twist, twin);
        }

        edges[edge].level = level;
        edges[edge].next = next;
        edges[edge].prev = prev;
        edges[edge].patch = patch;
      }

      edge_cnt++;
    }
    else if (tokens[0] == "P")
    {
      auto patch = patches.add(Patch({}));
      patches[patch].side = {(uint32_t)std::stoul(tokens[1]), 0};
    }
  }

  for (auto twin : parent_twins)
  {
    auto next = edges[edges[twin].next];
    auto next_parent_f = std::get_if<Parent>(&next.kind);
    if (next_parent_f)
    {
      auto edge_handles = edges[twin].handles();
      std::array<Id<Handle>, 2> twin_handles = {edge_handles.value()[1],
                                                edge_handles.value()[0]};
      Interpolant twist = edges[twin].twist;

      auto parent_twin = half_edge(next_parent_f->origin, twin_handles,
                                   Vector3(), twist, twin);
      edges[parent_twin].leftmost_child = edges[twin].twin;
      edges[twin].twin = parent_twin;
    }
  }

  for (auto stem : t_stems)
  {
    auto twin = edges[stem].twin;
    if (twin.has_value())
    {
      auto parent = edges[edges[edges[twin.value()].next].twin.value()].twin;

      edges[stem].kind =
          Child{parent.value(), {0.5, 0.5}, edges[stem].handles()};
    }
  }

  for (auto child : parentless)
  {
    auto parent = edges[edges[child].twin.value()].twin;
    if (parent.has_value())
    {
      Interval interval = {0.0, 0.5};

      edges[parent.value()].patch = edges[child].patch;

      auto child_next_twin = edges[edges[child].next].twin;
      if (child_next_twin.has_value())
      {
        auto next_child_edge = edges[edges[child_next_twin.value()].next].twin;
        if (next_child_edge.has_value() &&
            next_child_edge.value() != edges[child].twin.value())
        {
          edges[parent.value()].next = edges[child].next;
        }
      }
      else
      {
        interval = {0.5, 1.0};
        edges[parent.value()].next = edges[child].next;
      }

      auto child_prev_twin = edges[edges[child].prev].twin;
      if (child_prev_twin.has_value())
      {
        auto prev_child_edge = edges[edges[child_prev_twin.value()].prev].twin;
        if (prev_child_edge.has_value() &&
            prev_child_edge.value() != edges[child].twin.value())
        {
          edges[parent.value()].prev = edges[child].prev;
        }
      }
      else
      {
        interval = {0.0, 0.5};
        edges[parent.value()].prev = edges[child].prev;
      }

      edges[child].kind =
          Child{parent.value(), interval, edges[parent.value()].handles()};
    }
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
