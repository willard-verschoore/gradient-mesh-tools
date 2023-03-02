#include "gmt/half-edge.hpp"

#include "gmt/patch.hpp"

namespace gmt
{

Interval HalfEdge::interval() const
{
  return visit([](const Parent&) { return Interval(0.0f, 1.0f); },
               [](const Child& f) { return f.interval; }, *this);
}

Interval HalfEdge::relative_interval() const
{
  return is_parentable() ? Interval(0.0f, 1.0f) : interval();
}

bool HalfEdge::is_parentable() const
{
  return visit([=](const Parent&) { return true; },
               [=](const Child& f) { return is_empty(f.interval); }, *this);
}

std::optional<std::array<Id<Handle>, 2>> HalfEdge::handles() const
{
  return visit([](const Parent& a) { return std::make_optional(a.handles); },
               [](const Child& f) { return f.handles; }, *this);
}

void HalfEdge::set_handles(const std::array<Id<Handle>, 2>& handles)
{
  return visit([&](Parent& a) { a.handles = handles; },
               [&](Child& f) { f.handles = handles; }, *this);
}

void connect_twins(Id<HalfEdge> a, Id<HalfEdge> b, Storage<HalfEdge>& edges)
{
  edges[a].twin = b;
  edges[b].twin = a;
}

void connect(Id<HalfEdge> edge, Id<Patch> patch, Storage<HalfEdge>& edges,
             Storage<Patch>& patches)
{
  edges[edge].patch = patch;
  patches[patch].side = edge;
}

Id<HalfEdge> opposite(Id<HalfEdge> edge, const Storage<HalfEdge>& edges)
{
  return edges[edges[edge].next].next;
}

std::optional<Id<HalfEdge>> adjacent_next(Id<HalfEdge> edge,
                                          const Storage<HalfEdge>& edges)
{
  if (auto twin = edges[edges[edge].next].twin)
  {
    return edges[twin.value()].next;
  }
  else
  {
    return std::nullopt;
  }
}

std::optional<Id<HalfEdge>> adjacent_prev(Id<HalfEdge> edge,
                                          const Storage<HalfEdge>& edges)
{
  if (auto twin = edges[edges[edge].prev].twin)
  {
    return edges[twin.value()].prev;
  }
  else
  {
    return std::nullopt;
  }
}

Id<HalfEdge> eligible_parent(Id<HalfEdge> edge, const Storage<HalfEdge>& edges)
{
  return visit([=](const Parent&) { return edge; },
               [=](const Child& f)
               {
                 if (is_empty(f.interval))
                 {
                   return edge;
                 }
                 else
                 {
                   return f.parent;
                 }
               },
               edges[edge]);
}

ChildMap children(Id<HalfEdge> edge, const Storage<HalfEdge>& edges)
{
  ChildMap ret;
  auto child = edges[edge].leftmost_child;
  auto len = 1.0f;
  while (child.has_value() && len > 0.0f)
  {
    auto c = child.value();

    // Stop when the current child is not actually a child of the edge,
    // this can occur with t-junctions. Normally the length-part would
    // catch this, but floating point comparisons are finicky.
    auto child_f = std::get_if<Child>(&edges[c].kind);
    if (!child_f || child_f->parent != edge)
    {
      break;
    }

    auto i = edges[c].interval();
    ret.emplace(i.start, c);
    len -= length(i);
    child = adjacent_next(c, edges);
  }
  return ret;
}

ChildMap siblings(Id<HalfEdge> edge, const Storage<HalfEdge>& edges)
{
  auto parent = eligible_parent(edge, edges);
  if (parent != edge)
  {
    return children(parent, edges);
  }
  else
  {
    return {{0.0f, edge}};
  }
}

ChildMap abutting_edges(Id<HalfEdge> edge, const Storage<HalfEdge>& edges,
                        bool is_relative)
{
  ChildMap ret;
  const auto& e = edges[edge];
  auto interval = e.relative_interval();
  if (auto twin = e.twin)
  {
    auto twin_parent = eligible_parent(twin.value(), edges);
    for (auto [twin_start, twin_edge] : children(twin_parent, edges))
    {
      auto start = 1.0f - twin_start;
      if (contains(interval, start))
      {
        if (is_relative)
        {
          ret.emplace(relative_child_position(interval, start), twin_edge);
        }
        else
        {
          ret.emplace(start, twin_edge);
        }
      }
    }
  }
  return ret;
}

std::optional<Id<HalfEdge>> containing_twin(Id<HalfEdge> edge,
                                            const Storage<HalfEdge>& edges)
{
  const auto& e = edges[edge];
  auto interval = converse(e.relative_interval());
  if (auto twin = e.twin)
  {
    for (auto [twin_start, twin_edge] : siblings(twin.value(), edges))
    {
      auto twin_interval = edges[twin_edge].relative_interval();
      if (intersects(interval, twin_interval))
      {
        return twin_edge;
      }
    }
  }
  return std::nullopt;
}

std::size_t siblings_abutting_twin(Id<HalfEdge> edge,
                                   const Storage<HalfEdge>& edges)
{
  std::size_t total = 0;
  if (auto twin = containing_twin(edge, edges))
  {
    auto interval = converse(edges[twin.value()].relative_interval());
    for (auto [pos, sibling] : siblings(edge, edges))
    {
      auto child_interval = edges[sibling].relative_interval();
      if (intersects(interval, child_interval))
      {
        total++;
      }
    }
  }
  else
  {
    total++;
  }
  return total;
}

ChildIterator find(Id<HalfEdge> child, ChildMap& children)
{
  for (auto it = children.begin(); it != children.end(); ++it)
  {
    if (it->second == child)
    {
      return it;
    }
  }
  return children.end();
}

Id<HalfEdge> first_child(Id<HalfEdge> edge, const Storage<HalfEdge>& edges)
{
  if (edges[edge].leftmost_child.has_value())
  {
    return children(edge, edges).begin()->second;
  }
  return edge;
}

Id<HalfEdge> last_child(Id<HalfEdge> edge, const Storage<HalfEdge>& edges)
{
  if (edges[edge].leftmost_child.has_value())
  {
    return std::prev(children(edge, edges).end())->second;
  }
  return edge;
}

size_t num_children(Id<HalfEdge> edge, const Storage<HalfEdge>& edges)
{
  auto map = siblings(edge, edges);
  return map.size();
}

} // namespace gmt
