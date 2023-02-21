#include <QFile>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "gradient-mesh.hpp"

void read_cgm_handle(QStringList values, SlotMap<Handle>& handles)
{
  Id<HalfEdge> edge = {values[3].toUInt(), 0};
  Vector2 coord = Vector2(values[1].toFloat(), values[2].toFloat());

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

void read_cgm_point(QStringList values, std::vector<Id<HalfEdge>>& t_stems,
                    SlotMap<ControlPoint>& points)
{
  Id<HalfEdge> edge = {values[3].toUInt(), 0};
  Vector2 coord = Vector2(values[1].toFloat(), values[2].toFloat());

  auto point = points.add(ControlPoint({}));
  points[point].coords = coord;
  points[point].edge = edge;

  bool is_t_point = (values[5] == "1");
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

bool is_stem(Id<HalfEdge> edge, std::vector<Id<HalfEdge>>& t_stems)
{
  return std::find(t_stems.begin(), t_stems.end(), edge) != t_stems.end();
}

void GradientMesh::open_from_cgm(QString file_name)
{
  QFile newMesh(file_name);

  if (newMesh.open(QIODevice::ReadOnly))
  {
    QTextStream fileContents(&newMesh);

    QString currentLine;
    QStringList values;

    points.clear();
    handles.clear();
    edges.clear();
    patches.clear();

    std::set<Id<HalfEdge>> parent_twins;
    std::vector<Id<HalfEdge>> parentless;
    std::vector<Id<HalfEdge>> t_stems;

    uint edge_cnt = 0;

    while (!fileContents.atEnd())
    {
      currentLine = fileContents.readLine();
      values = currentLine.split(" ");

      if (values[0] == "S")
      {
        // Size of ControlPoints, HalfCurves and Patches
        // num_points = values[1].toUInt();
        // num_edges = values[2].toUInt();
        // num_patches = values[3].toUInt();
      }
      else if (values[0] == "H")
      {
        read_cgm_handle(values, handles);
      }
      else if (values[0] == "C")
      {
        read_cgm_point(values, t_stems, points);
      }
      else if (values[0] == "E")
      {
        Id<Patch> patch = {values[22].toUInt(), 0};

        // check if has twin we dont use boundary edges here
        if (patch.id < 65535)
        { // for some reason no patch is MAX_SHORT in pieters code
          Id<ControlPoint> corner = {values[1].toUInt(), 0};
          Id<HalfEdge> next = {values[19].toUInt(), 0};
          Id<HalfEdge> prev = {values[20].toUInt(), 0};

          std::optional<Id<HalfEdge>> twin = std::nullopt;
          if (values[21].toInt() >= 0)
          {
            twin = {values[21].toUInt(), 0};
          }

          Interpolant twist;
          twist.coords = Vector2(values[4].toFloat(), values[5].toFloat());
          twist.color = Vector3(values[15].toFloat(), values[16].toFloat(),
                                values[17].toFloat());

          Vector3 color = Vector3(values[6].toFloat(), values[7].toFloat(),
                                  values[8].toFloat());

          Id<Handle> handleTail = {values[2].toUInt(), 0};
          Id<Handle> handleHead = {values[3].toUInt(), 0};
          std::array<Id<Handle>, 2> edge_handles = {handleTail, handleHead};

          Id<HalfEdge> edge;
          int level = values[23].toInt();
          if (level > 0)
          {
            Interval interval = {0.0f, 1.0 / pow(2.0, double(level))};
            // Interval interval = {1.0 / pow(2.0, double(level)), 0.0};

            if (twin.has_value())
            {
              if (is_stem(twin.value(), t_stems))
              {
                edge = half_edge(corner, edge_handles, color, twist, twin);
                qDebug() << "twin is stem";
              }
              else
              {
                if (!is_stem({edge_cnt, 0}, t_stems))
                {
                  edge = half_edge(Id<HalfEdge>(), interval, edge_handles,
                                   color, twist, twin);

                  parentless.push_back(edge);
                  parent_twins.insert(twin.value());
                }
                else
                {
                  edge = half_edge(Id<HalfEdge>(), {0.5, 0.5}, edge_handles,
                                   color, twist, twin);
                  qDebug() << "created stem edge";
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
      else if (values[0] == "P")
      {
        auto patch = patches.add(Patch({}));
        patches[patch].side = {values[1].toUInt(), 0};
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

        qDebug() << "created parent edge";
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
      qDebug() << "connected stem";
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
          auto next_child_edge =
              edges[edges[child_next_twin.value()].next].twin;
          if (next_child_edge.has_value() &&
              next_child_edge.value() != edges[child].twin.value())
          {
            edges[parent.value()].next = edges[child].next;
            qDebug() << "set parent next";
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
          auto prev_child_edge =
              edges[edges[child_prev_twin.value()].prev].twin;
          if (prev_child_edge.has_value() &&
              prev_child_edge.value() != edges[child].twin.value())
          {
            edges[parent.value()].prev = edges[child].prev;
            qDebug() << "set parent prev";
          }
        }
        else
        {
          interval = {0.0, 0.5};
          edges[parent.value()].prev = edges[child].prev;
        }

        edges[child].kind =
            Child{parent.value(), interval, edges[parent.value()].handles()};
        qDebug() << "connected child to parent";
      }
      else
      {
        qDebug() << "could not determine parent for child edge";
      }
    }
  }
  qDebug() << "finished";
}

void GradientMesh::open_from_file(QString file_name)
{
  QFile file(file_name);

  points.clear();
  handles.clear();
  patches.clear();
  edges.clear();

  size_t num_points = 0;
  size_t num_handles = 0;
  size_t num_patches = 0;
  size_t num_edges = 0;

  if (file.open(QIODevice::ReadOnly))
  {
    QTextStream in(&file);
    QString line;
    QStringList tokens;

    line = in.readLine();
    line = in.readLine();
    tokens = line.split(" ", QString::SkipEmptyParts);
    num_points = tokens[0].toInt();
    num_handles = tokens[1].toInt();
    num_patches = tokens[2].toInt();
    num_edges = tokens[3].toInt();

    for (size_t i = 0; i < num_points; ++i)
    {
      line = in.readLine();
      if (line.isEmpty()) continue;
      tokens = line.split(" ", QString::SkipEmptyParts);

      auto point = points.add(ControlPoint({}));
      points[point].coords = Vector2(tokens[0].toFloat(), tokens[1].toFloat());
      points[point].edge = {tokens[2].toUInt(), 0};
    }

    for (size_t i = 0; i < num_handles; ++i)
    {
      line = in.readLine();
      if (line.isEmpty()) continue;
      tokens = line.split(" ", QString::SkipEmptyParts);

      auto handle = handles.add(Handle({}));
      handles[handle].edge = {tokens[0].toUInt(), 0};

      handles[handle].tangent.coords.x = tokens[1].toFloat();
      handles[handle].tangent.coords.y = tokens[2].toFloat();

      handles[handle].tangent.color.x = tokens[3].toFloat();
      handles[handle].tangent.color.y = tokens[4].toFloat();
      handles[handle].tangent.color.z = tokens[5].toFloat();
    }

    for (size_t i = 0; i < num_patches; ++i)
    {
      line = in.readLine();
      if (line.isEmpty()) continue;
      tokens = line.split(" ", QString::SkipEmptyParts);

      auto patch = patches.add(Patch({}));

      patches[patch].side = {tokens[0].toUInt(), 0};
    }

    for (size_t i = 0; i < num_edges; ++i)
    {
      line = in.readLine();
      if (line.isEmpty()) continue;
      tokens = line.split(" ", QString::SkipEmptyParts);

      Id<HalfEdge> edge;
      Interval interval = {tokens[0].toFloat(), tokens[1].toFloat()};
      Interpolant twist;
      twist.coords = Vector2(tokens[2].toFloat(), tokens[3].toFloat());
      twist.color = Vector3(tokens[4].toFloat(), tokens[5].toFloat(),
                            tokens[6].toFloat());
      Vector3 color = Vector3(tokens[7].toFloat(), tokens[8].toFloat(),
                              tokens[9].toFloat());

      std::optional<std::array<Id<Handle>, 2>> edge_handles = std::nullopt;
      if (tokens[10].toInt() >= 0)
      {
        Id<Handle> handle_1 = {tokens[10].toUInt(), 0};
        Id<Handle> handle_2 = {tokens[11].toUInt(), 0};
        edge_handles = {handle_1, handle_2};
      }

      std::optional<Id<HalfEdge>> twin = std::nullopt;
      if (tokens[12].toInt() >= 0)
      {
        twin = {tokens[12].toUInt(), 0};
      }
      Id<HalfEdge> prev = {tokens[13].toUInt(), 0};
      Id<HalfEdge> next = {tokens[14].toUInt(), 0};
      Id<Patch> patch = {tokens[15].toUInt(), 0};

      std::optional<Id<HalfEdge>> left_most_child = std::nullopt;
      if (tokens[16].toInt() >= 0)
      {
        left_most_child = {tokens[16].toUInt(), 0};
      }

      int is_child = tokens[17].toInt();
      if (is_child)
      {
        Id<HalfEdge> parent = {tokens[18].toUInt(), 0};
        std::optional<std::array<Id<Handle>, 2>> child_handles = std::nullopt;
        if (tokens[19].toInt() >= 0)
        {
          Id<Handle> handle_1 = {tokens[19].toUInt(), 0};
          Id<Handle> handle_2 = {tokens[20].toUInt(), 0};
          child_handles = {handle_1, handle_2};
        }

        float child_interval_start = tokens[21].toFloat();
        float child_interval_end = tokens[22].toFloat();
        Interval child_interval = {child_interval_start, child_interval_end};

        edge = half_edge(parent, child_interval, child_handles, color, twist,
                         twin);
      }
      else
      {
        Id<Handle> handle_1 = {tokens[18].toUInt(), 0};
        Id<Handle> handle_2 = {tokens[19].toUInt(), 0};
        Id<ControlPoint> origin = {tokens[20].toUInt(), 0};

        edge = half_edge(origin, {handle_1, handle_2}, color, twist, twin);
      }

      edges[edge].next = next;
      edges[edge].prev = prev;
      edges[edge].patch = patch;
      edges[edge].leftmost_child = left_most_child;
    }
  }

  file.close();
}

void GradientMesh::write_to_file(QTextStream& file) const
{
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
}
