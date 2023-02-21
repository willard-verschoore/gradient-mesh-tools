cmake_minimum_required(VERSION 3.10)
project(gradient-mesh-tools)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)

# Create a compile_commands.json file for language servers.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SOURCE_FILES
  src/draggable.cpp
  src/extrude.cpp
  src/gradient-mesh.cpp
  src/half-edge.cpp
  src/hermite.cpp
  src/interpolant.cpp
  src/io.cpp
  src/merge-conditions.cpp
  src/merging.cpp
  src/splitting.cpp
  src/tessellation.cpp
  src/vector2.cpp
  src/vector3.cpp
)

add_library(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC include)

find_package(Qt5 COMPONENTS Core Gui REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC Qt5::Core Qt5::Gui)

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
