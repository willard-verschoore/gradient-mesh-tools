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

void GradientMesh::open_from_file(std::string const &file_name)
{
  std::ifstream file(file_name);
  if (!file.is_open()) return;

  points.clear();
  handles.clear();
  patches.clear();
  edges.clear();

  size_t num_points = 0;
  size_t num_handles = 0;
  size_t num_patches = 0;
  size_t num_edges = 0;

  std::string line;
  std::vector<std::string> tokens;

  getline(file, line);
  getline(file, line);
  tokens = get_tokens(line, ' ');

  num_points = std::stoi(tokens[0]);
  num_handles = std::stoi(tokens[1]);
  num_patches = std::stoi(tokens[2]);
  num_edges = std::stoi(tokens[3]);

  for (size_t i = 0; i < num_points; ++i)
  {
    getline(file, line);
    if (line.empty()) continue;
    tokens = get_tokens(line, ' ');

    auto point = points.add(ControlPoint({}));
    points[point].coords = Vector2(std::stof(tokens[0]), std::stof(tokens[1]));
    points[point].edge = {(uint32_t)std::stoul(tokens[2]), 0};
  }

  for (size_t i = 0; i < num_handles; ++i)
  {
    getline(file, line);
    if (line.empty()) continue;
    tokens = get_tokens(line, ' ');

    auto handle = handles.add(Handle({}));
    handles[handle].edge = {(uint32_t)std::stoul(tokens[0]), 0};

    handles[handle].tangent.coords.x = std::stof(tokens[1]);
    handles[handle].tangent.coords.y = std::stof(tokens[2]);

    handles[handle].tangent.color.x = std::stof(tokens[3]);
    handles[handle].tangent.color.y = std::stof(tokens[4]);
    handles[handle].tangent.color.z = std::stof(tokens[5]);
  }

  for (size_t i = 0; i < num_patches; ++i)
  {
    getline(file, line);
    if (line.empty()) continue;
    tokens = get_tokens(line, ' ');

    auto patch = patches.add(Patch({}));
    patches[patch].side = {(uint32_t)std::stoul(tokens[0]), 0};
  }

  for (size_t i = 0; i < num_edges; ++i)
  {
    getline(file, line);
    if (line.empty()) continue;
    tokens = get_tokens(line, ' ');

    Id<HalfEdge> edge;
    Interval interval = {std::stof(tokens[0]), std::stof(tokens[1])};
    Interpolant twist;
    twist.coords = Vector2(std::stof(tokens[2]), std::stof(tokens[3]));
    twist.color = Vector3(std::stof(tokens[4]), std::stof(tokens[5]),
                          std::stof(tokens[6]));
    Vector3 color = Vector3(std::stof(tokens[7]), std::stof(tokens[8]),
                            std::stof(tokens[9]));

    std::optional<std::array<Id<Handle>, 2>> edge_handles = std::nullopt;
    if (std::stoi(tokens[10]) >= 0)
    {
      Id<Handle> handle_1 = {(uint32_t)std::stoul(tokens[10]), 0};
      Id<Handle> handle_2 = {(uint32_t)std::stoul(tokens[11]), 0};
      edge_handles = {handle_1, handle_2};
    }

    std::optional<Id<HalfEdge>> twin = std::nullopt;
    if (std::stoi(tokens[12]) >= 0)
    {
      twin = {(uint32_t)std::stoul(tokens[12]), 0};
    }
    Id<HalfEdge> prev = {(uint32_t)std::stoul(tokens[13]), 0};
    Id<HalfEdge> next = {(uint32_t)std::stoul(tokens[14]), 0};
    Id<Patch> patch = {(uint32_t)std::stoul(tokens[15]), 0};

    std::optional<Id<HalfEdge>> left_most_child = std::nullopt;
    if (std::stoi(tokens[16]) >= 0)
    {
      left_most_child = {(uint32_t)std::stoul(tokens[16]), 0};
    }

    int is_child = std::stoi(tokens[17]);
    if (is_child)
    {
      Id<HalfEdge> parent = {(uint32_t)std::stoul(tokens[18]), 0};
      std::optional<std::array<Id<Handle>, 2>> child_handles = std::nullopt;
      if (std::stoi(tokens[19]) >= 0)
      {
        Id<Handle> handle_1 = {(uint32_t)std::stoul(tokens[19]), 0};
        Id<Handle> handle_2 = {(uint32_t)std::stoul(tokens[20]), 0};
        child_handles = {handle_1, handle_2};
      }

      float child_interval_start = std::stof(tokens[21]);
      float child_interval_end = std::stof(tokens[22]);
      Interval child_interval = {child_interval_start, child_interval_end};

      edge =
          half_edge(parent, child_interval, child_handles, color, twist, twin);
    }
    else
    {
      Id<Handle> handle_1 = {(uint32_t)std::stoul(tokens[18]), 0};
      Id<Handle> handle_2 = {(uint32_t)std::stoul(tokens[19]), 0};
      Id<ControlPoint> origin = {(uint32_t)std::stoul(tokens[20]), 0};

      edge = half_edge(origin, {handle_1, handle_2}, color, twist, twin);
    }

    edges[edge].next = next;
    edges[edge].prev = prev;
    edges[edge].patch = patch;
    edges[edge].leftmost_child = left_most_child;
  }

  file.close();
}

void GradientMesh::write_to_file(std::string const &file_name) const
{
  std::ofstream file(file_name);
  if (!file.is_open()) return;

  file << "HEMESH\n";
  points.size();
  file << points.size() << " " << handles.size() << " " << patches.size() << " "
       << edges.size() << "\n";

  for (auto control_point : points)
  {
    file << control_point.coords.x << " ";
    file << control_point.coords.y << " ";
    file << control_point.edge.id << "\n";
  }

  for (auto handle : handles)
  {
    file << handle.edge.id << " ";
    file << handle.tangent.coords.x << " ";
    file << handle.tangent.coords.y << " ";

    file << handle.tangent.color.x << " ";
    file << handle.tangent.color.y << " ";
    file << handle.tangent.color.z << "\n";
  }

  for (auto patch : patches)
  {
    file << patch.side.id << "\n";
  }

  for (auto edge : edges)
  {
    file << edge.interval().start << " "; // 0
    file << edge.interval().end << " ";   // 1
    file << edge.twist.coords.x << " ";   // 2
    file << edge.twist.coords.y << " ";   // 3
    file << edge.twist.color.x << " ";    // 4
    file << edge.twist.color.y << " ";    // 5
    file << edge.twist.color.z << " ";    // 6
    file << edge.color.x << " ";          // 7
    file << edge.color.y << " ";          // 8
    file << edge.color.z << " ";          // 9

    if (edge.handles().has_value())
    {
      file << edge.handles().value()[0].id << " "; // 10
      file << edge.handles().value()[1].id << " "; // 1
    }
    else
    {
      file << -1 << " "; // 10
      file << -1 << " "; // 11
    }

    if (edge.twin.has_value())
    {
      file << edge.twin.value().id << " "; // 12
    }
    else
    {
      file << -1 << " "; // 12
    }
    file << edge.prev.id << " ";  // 13
    file << edge.next.id << " ";  // 14
    file << edge.patch.id << " "; // 15

    if (edge.leftmost_child.has_value())
    {
      file << edge.leftmost_child.value().id << " "; // 16
    }
    else
    {
      file << -1 << " "; // 16
    }

    auto child_f = std::get_if<Child>(&edge.kind);
    auto parent_f = std::get_if<Parent>(&edge.kind);
    if (child_f)
    {
      file << 1 << " ";                  // 17
      file << child_f->parent.id << " "; // 18
      if (child_f->handles.has_value())
      {
        file << child_f->handles.value()[0].id << " "; // 19
        file << child_f->handles.value()[1].id << " "; // 20
      }
      else
      {
        file << -1 << " "; // 19
        file << -1 << " "; // 20
      }
      file << child_f->interval.start << " "; // 21
      file << child_f->interval.end << " ";   // 22
    }
    else
    {
      file << 0 << " ";                       // 17
      file << parent_f->handles[0].id << " "; // 18
      file << parent_f->handles[1].id << " "; // 19
      file << parent_f->origin.id << " ";     // 20
    }

    file << "\n";
  }

  file.close();
}
