#include "gmt/gradient-mesh.hpp"

#include <catch2/catch.hpp>
#include <iostream>

using namespace gmt;
using namespace gmt::hermite;

static GradientMesh test_mesh()
{
  return GradientMesh{1.0f,
                      {Vector3(1.0f, 0.0f, 0.0f), Vector3(0.0f, 1.0f, 0.0f),
                       Vector3(0.0f, 0.0f, 1.0f), Vector3(1.0f, 0.0f, 1.0f)}};
}

TEST_CASE("Reading and writing patch data", "[hermite]")
{
  GradientMesh mesh = test_mesh();

  // Get patch data from mesh, alter it, and resubmit it.
  auto patch_data = mesh.patch_data();
  for (auto &patch : patch_data)
  {
    for (int i = 0; i < 16; ++i)
    {
      patch.data[i] = Interpolant{(float)std::rand() / (float)RAND_MAX,
                                  (float)std::rand() / (float)RAND_MAX,
                                  (float)std::rand() / (float)RAND_MAX,
                                  (float)std::rand() / (float)RAND_MAX,
                                  (float)std::rand() / (float)RAND_MAX};
    }
  }
  mesh.read_patch_data(patch_data, true);

  // Verify that the mesh received the submitted patch data.
  auto new_patch_data = mesh.patch_data();
  for (size_t i = 0; i < patch_data.size(); ++i)
    CHECK(patch_data[i] == new_patch_data[i]);
}
