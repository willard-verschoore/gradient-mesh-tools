cmake_minimum_required(VERSION 3.10)
project(gradient-mesh-tools)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)

# Create a compile_commands.json file for language servers.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SOURCE_FILES
  src/gmt/hermite/curve-matrix.cpp
  src/gmt/hermite/hermite.cpp
  src/gmt/hermite/interpolant.cpp
  src/gmt/hermite/patch-matrix.cpp
  src/gmt/hermite/vector2.cpp
  src/gmt/hermite/vector3.cpp
  src/gmt/extrude.cpp
  src/gmt/gradient-mesh.cpp
  src/gmt/half-edge.cpp
  src/gmt/io.cpp
  src/gmt/merge-conditions.cpp
  src/gmt/merging.cpp
  src/gmt/splitting.cpp
  src/gmt/tessellation.cpp
)

add_library(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC include)

# Set project to use standard (non-GNU) C++17.
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

# Enable warning flags, which differ between MSVC (Windows) and other compilers.
set(IS_CLANG "$<CXX_COMPILER_ID:Clang>")
set(IS_APPLE_CLANG "$<CXX_COMPILER_ID:AppleClang>")
set(IS_GNU "$<CXX_COMPILER_ID:GNU>")
set(IS_MSVC "$<CXX_COMPILER_ID:MSVC>")
target_compile_options(${PROJECT_NAME} PUBLIC
  $<$<OR:${IS_CLANG},${IS_APPLE_CLANG},${IS_GNU}>:
  -Wall
  -Wextra
  # Since we use mutually referencing data structures, we necessarily need to
  # set backreferences after construction.
  -Wno-missing-field-initializers
  >
  $<${IS_MSVC}:
  /W4
  >
)

# Enable unit tests and point CMake to the directory containing them.
enable_testing()
add_subdirectory(tests EXCLUDE_FROM_ALL)

# Generate documentation with Doxygen.
find_package(Doxygen)
if (DOXYGEN_FOUND)
  set(DOXYGEN_GENERATE_HTML YES)
  set(DOXYGEN_GENERATE_LATEX YES)
  set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/docs/")
  doxygen_add_docs(
    docs
    ${PROJECT_SOURCE_DIR}
    COMMENT "Generate Doxygen documentation"
  )
endif()
