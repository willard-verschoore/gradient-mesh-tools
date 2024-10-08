#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <iostream>

#include "gmt/gradient-mesh.hpp"

namespace gmt
{

using namespace hermite;

const int GradientMesh::MAX_SAMPLING_DENSITY = 64;

static bool equal_colors(Vector3 color1, Vector3 color2)
{
  Vector3 difference = 255.0f * (color1 - color2);
  return abs(difference.r) < 1.0f && abs(difference.g) < 1.0f &&
         abs(difference.b) < 1.0f;
}

hermite::Vector3 GradientMesh::find_color_opposite_origin(
    Id<HalfEdge> edge) const
{
  float t = edges[edge].interval().start;

  // If this edge is orthogonal to its parent as part of a t-junction we move up
  // to the parent edge before finding the twin.
  if (is_empty(edges[edge].interval()))
  {
    Child const *child = std::get_if<Child>(&edges[edge].kind);
    if (!child) return {-1.0f, -1.0f, -1.0f};
    edge = child->parent;
  }

  // Find the twin edge on the other side of the curve. Return an error value if
  // it does not exist.
  if (!edges[edge].twin.has_value()) return {-1.0f, -1.0f, -1.0f};
  Id<HalfEdge> opposite = edges[edge].twin.value();

  // Find the leaf node of the twin which overlaps with our edge's origin.
  Id<HalfEdge> old_opposite;
  while (old_opposite != opposite)
  {
    // Note that the opposing half edge goes in the opposite direction. So its
    // interval is inverted and we use 1 - t.
    old_opposite = opposite;
    opposite = find_child_contains_t(opposite, 1.0f - t);
  }

  CurveMatrix opposite_curve = curve_matrix(opposite);
  Interval opposite_interval = edges[opposite].relative_interval();
  float opposite_t = relative_child_position(opposite_interval, 1.0f - t);

  return interpolate(opposite_curve, opposite_t).color;
}

std::vector<Id<HalfEdge>> GradientMesh::find_continuous_t_junctions() const
{
  std::vector<Id<HalfEdge>> t_edges;

  for (auto it = edges.begin(); it != edges.end(); ++it)
  {
    // Every t-junction has an edge orthogonal to the head. This edge is a child
    // of the half edge along the bottom of the head and its interval is [t, t]
    // with t the coordinate of the t-junction along the parent.
    Child const *child = std::get_if<Child>(&it->kind);
    if (!child) continue;
    if (!is_empty(child->interval)) continue;

    Id<HalfEdge> edge = edges.get_handle(it);
    Vector3 opposite_color = find_color_opposite_origin(edge);

    // The t-junction is (partly) continuous if the color on one side of the
    // stem matches the color above the head. Alternatively, we could check
    // whether the origin color of the edge in the t-junction has not been
    // changed, i.e. it still equals the colour of its parent curve at t. There
    // are slight differences between these two approaches.
    if (equal_colors(it->color, opposite_color)) t_edges.push_back(edge);

    // We also check the color on the other side of the stem.
    if (!it->twin) continue;
    Id<HalfEdge> next_child = edges[it->twin.value()].next;
    if (equal_colors(edges[next_child].color, opposite_color))
      t_edges.push_back(next_child);
  }

  return t_edges;
}

void GradientMesh::fix_continuous_t_junctions(
    std::vector<Id<HalfEdge>> const &junctions)
{
  for (Id<HalfEdge> edge : junctions)
  {
    // An edge at a t-junction should always be a child.
    Child const *child = std::get_if<Child>(&edges[edge].kind);
    if (!child) continue;

    // Interpolate along the parent curve to get the new junction color.
    float t = child->interval.start;
    CurveMatrix parent_curve = curve_matrix(child->parent);
    Vector3 parent_color = interpolate(parent_curve, t).color;
    set_origin_color(edge, parent_color); // Traverses down the tree.
  }
}

std::vector<Vector3> GradientMesh::sample_colors(int density) const
{
  density = std::max(1, density);                    // Ensure positive density.
  density = std::min(density, MAX_SAMPLING_DENSITY); // Limit maximum density.

  std::vector<Vector3> samples;
  samples.reserve(patch_count() * (density + 1) * (density + 1));

  // Sample (density + 1)^2 colors from each patch.
  for (PatchMatrix const &patch : patch_data())
  {
    for (int x = 0; x <= density; ++x)
    {
      for (int y = 0; y <= density; ++y)
      {
        float u = (float)x / density;
        float v = (float)y / density;
        Interpolant sample = interpolate(patch, u, v);
        samples.push_back(sample.color);
      }
    }
  }

  return samples;
}

std::vector<float> GradientMesh::rgb_data(bool bezier, bool inactive) const
{
  std::vector<PatchMatrix> patches = patch_data();
  std::vector<float> data; // Stores resulting RGB data.
  int increment; // Matrix loop increment, 1 for all elements, 3 for corners.

  // All 16 elements of a Bezier patch matrix are valid RGB(XY) positions.
  // Each corresponds to a control point.
  if (bezier)
  {
    data.reserve(3 * 16 * patches.size()); // 3D RGB data, all 16 elements.
    increment = 1; // Used to loop over all matrix elements.
  }
  // Only the corners of a Hermite patch matrix are valid RGB(XY) positions.
  // They store the corners of the patch.
  else
  {
    data.reserve(3 * 4 * patches.size()); // 3D RGB data, only the 4 corners.
    increment = 3; // Used to loop over only corner matrix elements.
  }

  // Extract the RGB positions from all patches in the mesh.
  for (auto &patch : patches)
  {
    if (bezier) patch.to_bezier();

    for (int r = 0; r < 4; r += increment)
    {
      for (int c = 0; c < 4; c += increment)
      {
        data.push_back(patch(r, c).color.r);
        data.push_back(patch(r, c).color.g);
        data.push_back(patch(r, c).color.b);
      }
    }
  }

  if (!inactive) return data;

  // We can also include the data for "inactive" edges. These are edges which
  // have children running in parallel. The child edges contribute directly to a
  // patch and their data is therefore included in the previous step. The
  // underlying parent edge contributes indirectly to its children's patches. It
  // is therefore not included in the list of patch matrices, but it's data is
  // important because the children determine their origins and tangents based
  // on the parent's.
  for (auto it = edges.begin(); it != edges.end(); ++it)
  {
    Id<HalfEdge> edge = edges.get_handle(it);

    // Skip active, bottom level edges (leaves).
    if (children(edge, edges).empty()) continue;

    // A curve matrix looks like [origin, tangent, tangent, origin]. The
    // tangents are only valid RGB(XY) positions in Bezier form.
    CurveMatrix curve = curve_matrix(edge);
    if (bezier) curve.to_bezier();

    for (int i = 0; i < 4; i += increment)
    {
      data.push_back(curve[i].color.r);
      data.push_back(curve[i].color.g);
      data.push_back(curve[i].color.b);
    }
  }

  return data;
}

std::vector<float> GradientMesh::rgbxy_data(bool bezier, bool inactive) const
{
  std::vector<PatchMatrix> patches = patch_data();
  std::vector<float> data; // Stores resulting RGBXY data.
  int increment; // Matrix loop increment, 1 for all elements, 3 for corners.

  // All 16 elements of a Bezier patch matrix are valid RGBXY positions. Each
  // corresponds to a control point.
  if (bezier)
  {
    data.reserve(3 * 16 * patches.size()); // 5D RGBXY data, all 16 elements.
    increment = 1; // Used to loop over all matrix elements.
  }
  // Only the corners of a Hermite patch matrix are valid RGBXY positions. They
  // store the corners of the patch.
  else
  {
    data.reserve(3 * 4 * patches.size()); // 5D RGBXY data, only the 4 corners.
    increment = 3; // Used to loop over only corner matrix elements.
  }

  // Extract the RGBXY positions from all patches in the mesh.
  for (auto &patch : patches)
  {
    if (bezier) patch.to_bezier();

    for (int r = 0; r < 4; r += increment)
    {
      for (int c = 0; c < 4; c += increment)
      {
        data.push_back(patch(r, c).color.r);
        data.push_back(patch(r, c).color.g);
        data.push_back(patch(r, c).color.b);
        data.push_back(patch(r, c).coords.x);
        data.push_back(patch(r, c).coords.y);
      }
    }
  }

  if (!inactive) return data;

  // We can also include the data for "inactive" edges. These are edges which
  // have children running in parallel. The child edges contribute directly to a
  // patch and their data is therefore included in the previous step. The
  // underlying parent edge contributes indirectly to its children's patches. It
  // is therefore not included in the list of patch matrices, but it's data is
  // important because the children determine their origins and tangents based
  // on the parent's.
  for (auto it = edges.begin(); it != edges.end(); ++it)
  {
    Id<HalfEdge> edge = edges.get_handle(it);

    // Skip active, bottom level edges (leaves).
    if (children(edge, edges).empty()) continue;

    // A curve matrix looks like [origin, tangent, tangent, origin]. The
    // tangents are only valid RGBXY positions in Bezier form.
    CurveMatrix curve = curve_matrix(edge);
    if (bezier) curve.to_bezier();

    for (int i = 0; i < 4; i += increment)
    {
      data.push_back(curve[i].color.r);
      data.push_back(curve[i].color.g);
      data.push_back(curve[i].color.b);
      data.push_back(curve[i].coords.x);
      data.push_back(curve[i].coords.y);
    }
  }

  return data;
}

std::pair<std::vector<Vector3>, std::vector<uint32_t>>
GradientMesh::get_palette(size_t target_size, int sampling_density) const
{
  std::vector<Vector3> palette;
  std::vector<uint32_t> indices;

  if (!Py_IsInitialized())
  {
    Py_Initialize();
    _import_array(); // TODO: Check if < 0 like in import_array() macro.
  }

  // Load the Python module.
  std::string module_path = std::string(GMT_ROOT) + "/src";
  std::string module_load_command = "sys.path.append(\"" + module_path + "\")";
  PyRun_SimpleString("import sys");
  PyRun_SimpleString(module_load_command.c_str());

  PyObject *module_name, *module, *function, *arguments, *result;

  module_name = PyUnicode_FromString("decompose");
  module = PyImport_Import(module_name);
  Py_DECREF(module_name);
  if (module == NULL)
  {
    std::cout << "decompose can not be imported\n";
    PyErr_Print();
    return std::make_pair(palette, indices);
  }

  function = PyObject_GetAttrString(module, "get_palette");
  if (function == NULL || !PyCallable_Check(function))
  {
    std::cout << "get_palette is null or not callable\n";
    if (PyErr_Occurred()) PyErr_Print();
    Py_XDECREF(function);
    Py_DECREF(module);
    return std::make_pair(palette, indices);
  }

  std::vector<Vector3> rgb = sample_colors(sampling_density);
  npy_intp dims[2]{(npy_intp)rgb.size(), 3};
  PyArrayObject *np_rgb =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, dims, NPY_FLOAT, reinterpret_cast<void *>(rgb.data())));

  arguments = PyTuple_New(2);
  PyTuple_SetItem(arguments, 0, reinterpret_cast<PyObject *>(np_rgb));
  PyTuple_SetItem(arguments, 1, PyLong_FromSize_t(target_size));
  result = PyObject_CallObject(function, arguments);
  if (result == NULL)
  {
    std::cout << "function call went wrong\n";
    PyErr_Print();
    // Py_DECREF(result);
    // Py_DECREF(arguments);
    Py_DECREF(function);
    Py_DECREF(module);
    return std::make_pair(palette, indices);
  }

  PyArrayObject *np_palette =
      reinterpret_cast<PyArrayObject *>(PyTuple_GET_ITEM(result, 0));
  PyArrayObject *np_indices =
      reinterpret_cast<PyArrayObject *>(PyTuple_GET_ITEM(result, 1));

  Vector3 *raw_palette = reinterpret_cast<Vector3 *>(PyArray_DATA(np_palette));
  uint32_t *raw_indices =
      reinterpret_cast<uint32_t *>(PyArray_DATA(np_indices));

  uint32_t palette_size = PyArray_SIZE(np_palette) / 3;
  for (uint32_t i = 0; i < palette_size; ++i) palette.push_back(raw_palette[i]);

  uint32_t indices_size = PyArray_SIZE(np_indices);
  for (uint32_t i = 0; i < indices_size; ++i) indices.push_back(raw_indices[i]);

  Py_DECREF(function);
  Py_DECREF(module);
  Py_DECREF(np_rgb);
  Py_DECREF(np_palette);
  Py_DECREF(np_indices);
  // Py_FinalizeEx();

  return std::make_pair(palette, indices);
}

std::pair<std::vector<hermite::Vector3>, std::vector<uint32_t>>
GradientMesh::simplify_palette(
    std::vector<hermite::Vector3> const &palette) const
{
  std::vector<Vector3> simplified_palette;
  std::vector<uint32_t> simplified_indices;

  if (!Py_IsInitialized())
  {
    Py_Initialize();
    _import_array(); // TODO: Check if < 0 like in import_array() macro.
  }

  // Load the Python module.
  std::string module_path = std::string(GMT_ROOT) + "/src";
  std::string module_load_command = "sys.path.append(\"" + module_path + "\")";
  PyRun_SimpleString("import sys");
  PyRun_SimpleString(module_load_command.c_str());

  PyObject *module_name, *module, *function, *arguments, *result;

  module_name = PyUnicode_FromString("decompose");
  module = PyImport_Import(module_name);
  Py_DECREF(module_name);
  if (module == NULL)
  {
    std::cout << "decompose can not be imported\n";
    PyErr_Print();
    return std::make_pair(simplified_palette, simplified_indices);
  }

  function = PyObject_GetAttrString(module, "simplify_hull");
  if (function == NULL || !PyCallable_Check(function))
  {
    std::cout << "simplify_hull is null or not callable\n";
    if (PyErr_Occurred()) PyErr_Print();
    Py_XDECREF(function);
    Py_DECREF(module);
    return std::make_pair(simplified_palette, simplified_indices);
  }

  npy_intp palette_dims[2]{(npy_intp)palette.size(), 3};
  PyArrayObject *np_palette =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, palette_dims, NPY_FLOAT,
          reinterpret_cast<void *>(const_cast<Vector3 *>(palette.data()))));

  arguments = PyTuple_New(1);
  PyTuple_SetItem(arguments, 0, reinterpret_cast<PyObject *>(np_palette));
  result = PyObject_CallObject(function, arguments);
  if (result == NULL)
  {
    std::cout << "function call went wrong\n";
    PyErr_Print();
    // Py_DECREF(result);
    // Py_DECREF(arguments);
    Py_DECREF(function);
    Py_DECREF(module);
    return std::make_pair(simplified_palette, simplified_indices);
  }

  PyArrayObject *np_simplified_palette =
      reinterpret_cast<PyArrayObject *>(PyTuple_GET_ITEM(result, 0));
  PyArrayObject *np_simplified_indices =
      reinterpret_cast<PyArrayObject *>(PyTuple_GET_ITEM(result, 1));

  Vector3 *raw_simplified_palette =
      reinterpret_cast<Vector3 *>(PyArray_DATA(np_simplified_palette));
  uint32_t *raw_simplified_indices =
      reinterpret_cast<uint32_t *>(PyArray_DATA(np_simplified_indices));

  uint32_t simplified_palette_size = PyArray_SIZE(np_simplified_palette) / 3;
  for (uint32_t i = 0; i < simplified_palette_size; ++i)
    simplified_palette.push_back(raw_simplified_palette[i]);

  uint32_t simplified_indices_size = PyArray_SIZE(np_simplified_indices);
  for (uint32_t i = 0; i < simplified_indices_size; ++i)
    simplified_indices.push_back(raw_simplified_indices[i]);

  Py_DECREF(function);
  Py_DECREF(module);
  Py_DECREF(np_simplified_palette);
  Py_DECREF(np_simplified_indices);
  // Py_FinalizeEx();

  return std::make_pair(simplified_palette, simplified_indices);
}

std::vector<hermite::Vector3> GradientMesh::optimize_palette(
    std::vector<hermite::Vector3> const &palette, double lambda,
    int sampling_density) const
{
  std::vector<Vector3> optimized_palette;

  if (!Py_IsInitialized())
  {
    Py_Initialize();
    _import_array(); // TODO: Check if < 0 like in import_array() macro.
  }

  // Load the Python module.
  std::string module_path = std::string(GMT_ROOT) + "/src";
  std::string module_load_command = "sys.path.append(\"" + module_path + "\")";
  PyRun_SimpleString("import sys");
  PyRun_SimpleString(module_load_command.c_str());

  PyObject *module_name, *module, *function, *arguments, *result;

  module_name = PyUnicode_FromString("decompose");
  module = PyImport_Import(module_name);
  Py_DECREF(module_name);
  if (module == NULL)
  {
    std::cout << "decompose can not be imported\n";
    PyErr_Print();
    return optimized_palette;
  }

  function = PyObject_GetAttrString(module, "optimize_vertex_positions");
  if (function == NULL || !PyCallable_Check(function))
  {
    std::cout << "optimize_vertex_positions is null or not callable\n";
    if (PyErr_Occurred()) PyErr_Print();
    Py_XDECREF(function);
    Py_DECREF(module);
    return optimized_palette;
  }

  npy_intp palette_dims[2]{(npy_intp)palette.size(), 3};
  PyArrayObject *np_palette =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, palette_dims, NPY_FLOAT,
          reinterpret_cast<void *>(const_cast<Vector3 *>(palette.data()))));

  std::vector<Vector3> rgb = sample_colors(sampling_density);
  npy_intp rgb_dims[2]{(npy_intp)rgb.size(), 3};
  PyArrayObject *np_rgb =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, rgb_dims, NPY_FLOAT, reinterpret_cast<void *>(rgb.data())));

  arguments = PyTuple_New(3);
  PyTuple_SetItem(arguments, 0, reinterpret_cast<PyObject *>(np_palette));
  PyTuple_SetItem(arguments, 1, reinterpret_cast<PyObject *>(np_rgb));
  PyTuple_SetItem(arguments, 2, PyFloat_FromDouble(lambda));
  result = PyObject_CallObject(function, arguments);
  if (result == NULL)
  {
    std::cout << "function call went wrong\n";
    PyErr_Print();
    // Py_DECREF(result);
    // Py_DECREF(arguments);
    Py_DECREF(function);
    Py_DECREF(module);
    return optimized_palette;
  }

  PyArrayObject *np_optimized_palette =
      reinterpret_cast<PyArrayObject *>(result);

  Vector3 *raw_optimized_palette =
      reinterpret_cast<Vector3 *>(PyArray_DATA(np_optimized_palette));

  uint32_t optimized_palette_size = PyArray_SIZE(np_optimized_palette) / 3;
  for (uint32_t i = 0; i < optimized_palette_size; ++i)
    optimized_palette.push_back(raw_optimized_palette[i]);

  Py_DECREF(function);
  Py_DECREF(module);
  Py_DECREF(np_rgb);
  Py_DECREF(np_optimized_palette);
  // Py_FinalizeEx();

  return optimized_palette;
}

std::vector<float> GradientMesh::get_weights(
    std::vector<Vector3> const &palette, WeightType type, bool bezier,
    bool inactive) const
{
  if (!Py_IsInitialized())
  {
    Py_Initialize();
    _import_array(); // TODO: Check if < 0 like in import_array() macro.
  }

  // Load the Python module.
  std::string module_path = std::string(GMT_ROOT) + "/src";
  std::string module_load_command = "sys.path.append(\"" + module_path + "\")";
  PyRun_SimpleString("import sys");
  PyRun_SimpleString(module_load_command.c_str());

  PyObject *module_name, *module, *function, *arguments, *result;

  module_name = PyUnicode_FromString("decompose");
  module = PyImport_Import(module_name);
  Py_DECREF(module_name);
  if (module == NULL)
  {
    std::cout << "decompose can not be imported\n";
    PyErr_Print();
    return std::vector<float>(0);
  }

  function = PyObject_GetAttrString(module, "get_weights");
  if (function == NULL || !PyCallable_Check(function))
  {
    std::cout << "get_weights is null or not callable\n";
    if (PyErr_Occurred()) PyErr_Print();
    Py_XDECREF(function);
    Py_DECREF(module);
    return std::vector<float>(0);
  }

  std::vector<float> rgbxy = rgbxy_data(bezier, inactive);
  npy_intp dims[2]{(npy_intp)rgbxy.size() / 5, 5};
  PyArrayObject *np_rgbxy =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, dims, NPY_FLOAT, reinterpret_cast<void *>(rgbxy.data())));

  npy_intp p_dims[2]{(npy_intp)palette.size(), 3};
  PyArrayObject *np_palette =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, p_dims, NPY_FLOAT,
          reinterpret_cast<void *>(const_cast<Vector3 *>(palette.data()))));

  arguments = PyTuple_New(3);
  PyTuple_SetItem(arguments, 0, reinterpret_cast<PyObject *>(np_rgbxy));
  PyTuple_SetItem(arguments, 1, reinterpret_cast<PyObject *>(np_palette));
  PyTuple_SetItem(arguments, 2, PyLong_FromLong((long)type));
  result = PyObject_CallObject(function, arguments);
  if (result == NULL)
  {
    std::cout << "function call went wrong\n";
    PyErr_Print();
    // Py_DECREF(result);
    // Py_DECREF(arguments);
    Py_DECREF(function);
    Py_DECREF(module);
    return std::vector<float>(0);
  }

  PyArrayObject *np_weights = reinterpret_cast<PyArrayObject *>(result);
  float *raw_weights = reinterpret_cast<float *>(PyArray_DATA(np_weights));
  int weights_size = PyArray_SIZE(np_weights);

  std::vector<float> weights(weights_size);
  for (int i = 0; i < weights_size; ++i) weights[i] = raw_weights[i];

  Py_DECREF(function);
  Py_DECREF(module);
  Py_DECREF(np_rgbxy);
  Py_DECREF(np_weights);
  // Py_FinalizeEx();

  return weights;
}

void GradientMesh::recolor(std::vector<float> const &weights,
                           std::vector<Vector3> const &palette, bool bezier,
                           bool inactive, bool create_tangents)
{
  std::vector<Id<HalfEdge>> junctions = find_continuous_t_junctions();
  std::vector<PatchMatrix> patches = patch_data();
  int increment = bezier ? 1 : 3;     // 1 for all elements, 3 for just corners.
  Vector3 color = {0.0f, 0.0f, 0.0f}; // Stores a sum of palette colors.
  size_t index = 0;                   // Tracks the current weight vector.

  // Recolor the all patches in the mesh. In the Bezier case we recolor all
  // matrix elements while in the Hermite case we only recolor the corners.
  for (auto &patch : patches)
  {
    // Convert Hermite patches to Bezier for recoloring if needed.
    if (bezier) patch.to_bezier();

    for (int r = 0; r < 4; r += increment)
    {
      for (int c = 0; c < 4; c += increment)
      {
        color = {0.0f, 0.0f, 0.0f};

        // Each color is a weighted sum of palette colors.
        for (size_t p = 0; p < palette.size(); ++p)
          color += weights[index * palette.size() + p] * palette[p];
        patch(r, c).color = color;

        ++index;
      }
    }

    // Convert Bezier patches back to Hermite for submission to the mesh.
    if (bezier) patch.to_hermite();
  }

  // Send the recolored patches to the mesh.
  // TODO: Is this function necessary? Use read_patch_matrix in the loop above?
  read_patch_data(patches, create_tangents);

  if (!inactive)
  {
    fix_inactive_origin_colors();
    fix_continuous_t_junctions(junctions);
    return;
  }

  // Recolor inactive edges (edges with parallel children).
  for (auto it = edges.begin(); it != edges.end(); ++it)
  {
    Id<HalfEdge> edge = edges.get_handle(it);

    // Skip active, bottom level edges (leaves).
    if (children(edge, edges).empty()) continue;

    CurveMatrix curve = curve_matrix(edge);

    // Convert Hermite curves to Bezier for recoloring if needed.
    if (bezier) curve.to_bezier();

    for (int i = 0; i < 4; i += increment)
    {
      color = {0.0f, 0.0f, 0.0f};

      // Each color is a weighted sum of palette colors.
      for (size_t p = 0; p < palette.size(); ++p)
        color += weights[index * palette.size() + p] * palette[p];
      curve[i].color = color;

      ++index;
    }

    // Convert Bezier curves back to Hermite for submission to the mesh.
    if (bezier) curve.to_hermite();

    read_curve_matrix(edge, curve, create_tangents);
  }

  fix_inactive_origin_colors();
  fix_continuous_t_junctions(junctions);
}

} // namespace gmt
