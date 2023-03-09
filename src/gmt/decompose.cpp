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
  std::vector<float> data;
  data.reserve(3 * edges.size());

  for (auto const &edge : edges)
  {
    data.push_back(edge.color.r);
    data.push_back(edge.color.g);
    data.push_back(edge.color.b);
  }

  return data;
}

std::vector<float> GradientMesh::rgbxy_data() const
{
  std::vector<float> data;
  data.reserve(5 * edges.size());

  for (auto const &edge : edges)
  {
    data.push_back(edge.color.r);
    data.push_back(edge.color.g);
    data.push_back(edge.color.b);

    // Find the origin position of the edge. Easy for parent edges, but
    // more difficult for child edges.
    Vector2 origin;
    visit([&](Parent const &parent) { origin = points[parent.origin].coords; },
          [&](Child const &child)
          {
            CurveMatrix curve = curve_matrix(child.parent);
            float t = child.interval.start;
            origin = interpolate(curve, t).coords;
          },
          edge);

    data.push_back(origin.x);
    data.push_back(origin.y);
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
  npy_intp dims[2]{edges.size(), 3};
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
  npy_intp dims[2]{edges.size(), 5};
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
  int palette_size = palette.size();
  int i = 0;
  for (auto &edge : edges)
  {
    // Each color is a weighted sum of palette colors.
    edge.color = {0.0f, 0.0f, 0.0f};
    for (int j = 0; j < palette_size; ++j)
      edge.color += weights[i * palette_size + j] * palette[j];
    ++i;
  }
}

} // namespace gmt
