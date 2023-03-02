#include "gmt/gradient-mesh.hpp"

// This file holds functions related to the merge conditions.

namespace gmt
{

/// Checks for c0-continuity at the edge between left and right.
inline bool is_c0(hermite::PatchMatrix left, hermite::PatchMatrix right,
                  MergeOptions merge_options)
{
  auto r1 = (left(0, 3)) - (right(0, 0));
  auto r2 = (left(0, 3) + left(1, 3) / 3) - (right(0, 0) + right(1, 0) / 3);
  auto r3 = (left(3, 3) - left(2, 3) / 3) - (right(3, 0) - right(2, 0) / 3);
  auto r4 = (left(3, 3)) - (right(3, 0));

  return is_zero(r1, merge_options.coords.c0, merge_options.color.c0) &&
         is_zero(r2, merge_options.coords.c0, merge_options.color.c0) &&
         is_zero(r3, merge_options.coords.c0, merge_options.color.c0) &&
         is_zero(r4, merge_options.coords.c0, merge_options.color.c0);
}

/// Checks for c1-continuity (and thus also c0) at the edge between left and
/// right.
inline bool is_c1(hermite::PatchMatrix left, hermite::PatchMatrix right,
                  float t, MergeOptions merge_options)
{
  if (!is_c0(left, right, merge_options))
  {
    return false;
  }

  auto t_left = 1.0f - t;
  auto t_right = t;
  auto r1 = (left(0, 2) * t_left) - (right(0, 1) * t_right);
  auto r2 = ((left(0, 2) + left(1, 2) / 3) * t_left) -
            ((right(0, 1) + right(1, 1) / 3) * t_right);
  auto r3 = ((left(3, 2) - left(2, 2) / 3) * t_left) -
            ((right(3, 1) - right(2, 1) / 3) * t_right);
  auto r4 = (left(3, 2) * t_left) - (right(3, 1) * t_right);

  return is_zero(r1, merge_options.coords.c1, merge_options.color.c1) &&
         is_zero(r2, merge_options.coords.c1, merge_options.color.c1) &&
         is_zero(r3, merge_options.coords.c1, merge_options.color.c1) &&
         is_zero(r4, merge_options.coords.c1, merge_options.color.c1);
}

/// Checks for c2-continuity (and thus also c1, c0) at the edge between left and
/// right.
inline bool is_c2(hermite::PatchMatrix left, hermite::PatchMatrix right,
                  float t, MergeOptions merge_options)
{
  if (!is_c1(left, right, t, merge_options))
  {
    return false;
  }

  auto t_left = pow(1.0f - t, 2);
  auto t_right = pow(t, 2);
  auto r1 =
      ((6 * left(0, 0) + 2 * left(0, 1) + 4 * left(0, 2) - 6 * left(0, 3)) *
       t_left) -
      ((-6 * right(0, 0) - 4 * right(0, 1) - 2 * right(0, 2) +
        6 * right(0, 3)) *
       t_right);
  auto r2 =
      ((6 * left(0, 0) + 2 * left(1, 0) + 2 * left(0, 1) + 2 * left(1, 1) / 3 +
        4 * left(0, 2) + 4 * left(1, 2) / 3 - 6 * left(0, 3) - 2 * left(1, 3)) *
       t_left) -
      ((-6 * right(0, 0) - 2 * right(1, 0) - 4 * right(0, 1) -
        4 * right(1, 1) / 3 - 2 * right(0, 2) - 2 * right(1, 2) / 3 +
        6 * right(0, 3) + 2 * right(1, 3)) *
       t_right);
  auto r3 =
      ((-2 * left(2, 0) + 6 * left(3, 0) - 2 * left(2, 1) / 3 + 2 * left(3, 1) -
        4 * left(2, 2) / 3 + 4 * left(3, 2) + 2 * left(2, 3) - 6 * left(3, 3)) *
       t_left) -
      ((2 * right(2, 0) - 6 * right(3, 0) + 4 * right(2, 1) / 3 -
        4 * right(3, 1) + 2 * right(2, 2) / 3 - 2 * right(3, 2) -
        2 * right(2, 3) + 6 * right(3, 3)) *
       t_right);
  auto r4 =
      ((6 * left(3, 0) + 2 * left(3, 1) + 4 * left(3, 2) - 6 * left(3, 3)) *
       t_left) -
      ((-6 * right(3, 0) - 4 * right(3, 1) - 2 * right(3, 2) +
        6 * right(3, 3)) *
       t_right);

  return is_zero(r1, merge_options.coords.c2, merge_options.color.c2) &&
         is_zero(r2, merge_options.coords.c2, merge_options.color.c2) &&
         is_zero(r3, merge_options.coords.c2, merge_options.color.c2) &&
         is_zero(r4, merge_options.coords.c2, merge_options.color.c2);
}

/// Checks for c3-continuity (and thus also c2, c1, c0) at the edge between left
/// and right.
inline bool is_c3(hermite::PatchMatrix left, hermite::PatchMatrix right,
                  float t, MergeOptions merge_options)
{
  if (!is_c2(left, right, t, merge_options))
  {
    return false;
  }

  auto t_left = pow(1.0f - t, 3);
  auto t_right = pow(t, 3);
  auto r1 =
      ((12 * left(0, 0) + 6 * left(0, 1) + 6 * left(0, 2) - 12 * left(0, 3)) *
       t_left) -
      ((12 * right(0, 0) + 6 * right(0, 1) + 6 * right(0, 2) -
        12 * right(0, 3)) *
       t_right);
  auto r2 =
      ((12 * left(0, 0) + 4 * left(1, 0) + 6 * left(0, 1) + 2 * left(1, 1) +
        6 * left(0, 2) + 2 * left(1, 2) - 12 * left(0, 3) - 4 * left(1, 3)) *
       t_left) -
      ((12 * right(0, 0) + 4 * right(1, 0) + 6 * right(0, 1) + 2 * right(1, 1) +
        6 * right(0, 2) + 2 * right(1, 2) - 12 * right(0, 3) -
        4 * right(1, 3)) *
       t_right);
  auto r3 =
      ((-4 * left(2, 0) + 12 * left(3, 0) - 2 * left(2, 1) + 6 * left(3, 1) -
        2 * left(2, 2) + 6 * left(3, 2) + 4 * left(2, 3) - 12 * left(3, 3)) *
       t_left) -
      ((-4 * right(2, 0) + 12 * right(3, 0) - 2 * right(2, 1) +
        6 * right(3, 1) - 2 * right(2, 2) + 6 * right(3, 2) + 4 * right(2, 3) -
        12 * right(3, 3)) *
       t_right);
  auto r4 =
      ((12 * left(3, 0) + 6 * left(3, 1) + 6 * left(3, 2) - 12 * left(3, 3)) *
       t_left) -
      ((12 * right(3, 0) + 6 * right(3, 1) + 6 * right(3, 2) -
        12 * right(3, 3)) *
       t_right);

  return is_zero(r1, merge_options.coords.c3, merge_options.color.c3) &&
         is_zero(r2, merge_options.coords.c3, merge_options.color.c3) &&
         is_zero(r3, merge_options.coords.c3, merge_options.color.c3) &&
         is_zero(r4, merge_options.coords.c3, merge_options.color.c3);
}

float GradientMesh::potential_splitting_factor(hermite::PatchMatrix left,
                                               hermite::PatchMatrix right)
{
  return length(left(0, 2) / 3) / length((right(0, 1) + left(0, 2)) / 3);
}

bool GradientMesh::can_merge_neighbours(Id<HalfEdge> separator_edge,
                                        MergeOptions merge_options)
{
  if (!edges[separator_edge].twin.has_value())
  {
    return false;
  }

  auto top_left = edges[separator_edge].prev;
  auto potential_top_right = adjacent_next(top_left, edges);
  assert(potential_top_right.has_value());
  auto top_right = potential_top_right.value();

  // Since the top_right is determined using adjacent_next, it can be a parent.
  // In that case, we should use the child, otherwise determining the patch
  // matrix is very difficult.
  top_right = first_child(top_right, edges);

  auto hermite_left = patch_matrix(top_left);
  auto hermite_right = patch_matrix(top_right);
  auto t = potential_splitting_factor(hermite_left, hermite_right);

  return is_c3(hermite_left, hermite_right, t, merge_options);
}

} // namespace gmt
