find_package(Catch2)

# Only build tests if the Catch library is installed on the system.
if (Catch2_FOUND)
  set(SOURCE_FILES
    catch.cpp
    gradient-mesh.cpp
    hermite.cpp
    slot-map.cpp
  )

  add_executable(tests ${SOURCE_FILES})
  target_link_libraries(tests PRIVATE ${PROJECT_NAME} Catch2::Catch2)

  include(CTest)
  include(Catch)
  catch_discover_tests(tests)

else ()
  message(WARNING "The Catch2 library was not found on the system, so package tests will not be built.")
endif(Catch2_FOUND)
