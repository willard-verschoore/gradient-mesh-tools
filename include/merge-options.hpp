#pragma once

/**
 * These structs specify the allowed errors in the continuity conditions,
 * split in coordinates & color and by continuity level.
 */

/// Allowed errors for c0, c1, c2, and c3 continuity during merging.
struct ContinuityErrors
{
  float c0, c1, c2, c3;
};

/// Specifies allowed continuity errors for the coordinates and colors
/// respectively.
struct MergeOptions
{
  ContinuityErrors coords = {0.00001f, 0.00001f, 0.00001f, 0.00001f};
  ContinuityErrors color = {0.00001f, 0.00001f, 0.00001f, 0.00001f};
};
