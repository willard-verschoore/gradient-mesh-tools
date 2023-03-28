#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <iostream>

#include "gmt/gradient-mesh.hpp"

namespace gmt
{

using namespace hermite;

std::vector<float> GradientMesh::rgb_data() const
{
  std::vector<PatchMatrix> patches = patch_data();
  std::vector<float> data;
  data.reserve(3 * 4 * patches.size()); // 3D RGB data, 4 corners.

  for (auto const &patch : patches)
  {
    // The corners of a patch matrix store the corners of the patch. The order
    // of the corners is chosen such that it matches the recolor() function.
    data.push_back(patch(0, 0).color.r);
    data.push_back(patch(0, 0).color.g);
    data.push_back(patch(0, 0).color.b);
    data.push_back(patch(3, 0).color.r);
    data.push_back(patch(3, 0).color.g);
    data.push_back(patch(3, 0).color.b);
    data.push_back(patch(3, 3).color.r);
    data.push_back(patch(3, 3).color.g);
    data.push_back(patch(3, 3).color.b);
    data.push_back(patch(0, 3).color.r);
    data.push_back(patch(0, 3).color.g);
    data.push_back(patch(0, 3).color.b);
  }

  return data;
}

std::vector<float> GradientMesh::rgbxy_data() const
{
  std::vector<PatchMatrix> patches = patch_data();
  std::vector<float> data;
  data.reserve(5 * 4 * patches.size()); // 5D RGBXY data, 4 corners.

  for (auto const &patch : patches)
  {
    // The corners of a patch matrix store the corners of the patch. The order
    // of the corners is chosen such that it matches the recolor() function.
    data.push_back(patch(0, 0).color.r);
    data.push_back(patch(0, 0).color.g);
    data.push_back(patch(0, 0).color.b);
    data.push_back(patch(0, 0).coords.x);
    data.push_back(patch(0, 0).coords.y);
    data.push_back(patch(3, 0).color.r);
    data.push_back(patch(3, 0).color.g);
    data.push_back(patch(3, 0).color.b);
    data.push_back(patch(3, 0).coords.x);
    data.push_back(patch(3, 0).coords.y);
    data.push_back(patch(3, 3).color.r);
    data.push_back(patch(3, 3).color.g);
    data.push_back(patch(3, 3).color.b);
    data.push_back(patch(3, 3).coords.x);
    data.push_back(patch(3, 3).coords.y);
    data.push_back(patch(0, 3).color.r);
    data.push_back(patch(0, 3).color.g);
    data.push_back(patch(0, 3).color.b);
    data.push_back(patch(0, 3).coords.x);
    data.push_back(patch(0, 3).coords.y);
  }

  return data;
}

std::vector<Vector3> GradientMesh::get_palette() const
{
  std::vector<Vector3> palette;

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
    return palette;
  }

  function = PyObject_GetAttrString(module, "get_palette");
  if (function == NULL || !PyCallable_Check(function))
  {
    std::cout << "get_palette is null or not callable\n";
    if (PyErr_Occurred()) PyErr_Print();
    Py_XDECREF(function);
    Py_DECREF(module);
    return palette;
  }

  std::vector<float> rgb = rgb_data();
  npy_intp dims[2]{(npy_intp)rgb.size() / 3, 3};
  PyArrayObject *np_rgb =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, dims, NPY_FLOAT, reinterpret_cast<void *>(rgb.data())));

  arguments = PyTuple_New(1);
  PyTuple_SetItem(arguments, 0, reinterpret_cast<PyObject *>(np_rgb));
  result = PyObject_CallObject(function, arguments);
  if (result == NULL)
  {
    std::cout << "function call went wrong\n";
    PyErr_Print();
    // Py_DECREF(result);
    // Py_DECREF(arguments);
    Py_DECREF(function);
    Py_DECREF(module);
    return palette;
  }

  PyArrayObject *np_palette = reinterpret_cast<PyArrayObject *>(result);
  Vector3 *raw_palette = reinterpret_cast<Vector3 *>(PyArray_DATA(np_palette));
  int palette_size = PyArray_SIZE(np_palette) / 3;
  for (int i = 0; i < palette_size; ++i) palette.push_back(raw_palette[i]);

  Py_DECREF(function);
  Py_DECREF(module);
  Py_DECREF(np_rgb);
  Py_DECREF(np_palette);
  // Py_FinalizeEx();

  return palette;
}

std::vector<uint32_t> GradientMesh::get_palette_indices() const
{
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
    return indices;
  }

  function = PyObject_GetAttrString(module, "get_palette_indices");
  if (function == NULL || !PyCallable_Check(function))
  {
    std::cout << "get_palette_indices is null or not callable\n";
    if (PyErr_Occurred()) PyErr_Print();
    Py_XDECREF(function);
    Py_DECREF(module);
    return indices;
  }

  std::vector<float> rgb = rgb_data();
  npy_intp dims[2]{(npy_intp)rgb.size() / 3, 3};
  PyArrayObject *np_rgb =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, dims, NPY_FLOAT, reinterpret_cast<void *>(rgb.data())));

  arguments = PyTuple_New(1);
  PyTuple_SetItem(arguments, 0, reinterpret_cast<PyObject *>(np_rgb));
  result = PyObject_CallObject(function, arguments);
  if (result == NULL)
  {
    std::cout << "function call went wrong\n";
    PyErr_Print();
    // Py_DECREF(result);
    // Py_DECREF(arguments);
    Py_DECREF(function);
    Py_DECREF(module);
    return indices;
  }

  PyArrayObject *np_indices = reinterpret_cast<PyArrayObject *>(result);
  uint32_t *raw_indices =
      reinterpret_cast<uint32_t *>(PyArray_DATA(np_indices));
  uint32_t indices_size = PyArray_SIZE(np_indices);
  for (uint32_t i = 0; i < indices_size; ++i) indices.push_back(raw_indices[i]);

  Py_DECREF(function);
  Py_DECREF(module);
  Py_DECREF(np_rgb);
  Py_DECREF(np_indices);
  // Py_FinalizeEx();

  return indices;
}

std::vector<float> GradientMesh::get_weights(
    std::vector<Vector3> const &palette) const
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

  std::vector<float> rgbxy = rgbxy_data();
  npy_intp dims[2]{(npy_intp)rgbxy.size() / 5, 5};
  PyArrayObject *np_rgbxy =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, dims, NPY_FLOAT, reinterpret_cast<void *>(rgbxy.data())));

  npy_intp p_dims[2]{(npy_intp)palette.size(), 3};
  PyArrayObject *np_palette =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNewFromData(
          2, p_dims, NPY_FLOAT,
          reinterpret_cast<void *>(const_cast<Vector3 *>(palette.data()))));

  arguments = PyTuple_New(2);
  PyTuple_SetItem(arguments, 0, reinterpret_cast<PyObject *>(np_rgbxy));
  PyTuple_SetItem(arguments, 1, reinterpret_cast<PyObject *>(np_palette));
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
                           std::vector<Vector3> const &palette)
{
  size_t index = 0;
  Vector3 color = {0.0f, 0.0f, 0.0f};

  for (auto const &patch : patches)
  {
    // Get the sides of the patch.
    Id<HalfEdge> sides[4];
    sides[0] = patch.side;
    sides[1] = edges[sides[0]].next;
    sides[2] = edges[sides[1]].next;
    sides[3] = edges[sides[2]].next;

    for (size_t i = 0; i < 4; ++i)
    {
      // Each color is a weighted sum of palette colors.
      color = {0.0f, 0.0f, 0.0f};
      for (size_t j = 0; j < palette.size(); ++j)
        color += weights[index * palette.size() + j] * palette[j];
      set_color(sides[i], color);
      ++index;
    }
  }
}

} // namespace gmt
