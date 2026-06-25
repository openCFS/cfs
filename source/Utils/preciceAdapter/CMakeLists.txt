# Define the preciceadapter library
add_library(preciceadapter STATIC)

# Set the library's include directories to the current source directory
target_include_directories(preciceadapter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


if(USE_PRECICE)
  # Add the real adapter sources (independent of how preCICE is provided)
  target_sources(preciceadapter PRIVATE
  MinimalXmlParser.cc
  PreciceAdapter.cc
  IPreciceAdapter.cc
  NodeResult.cc
  ElementResult.cc
  CoefFunctionGridElemDefaultPrecice.cc
  CoefFunctionGridNodalDefaultPrecice.cc
  TransientDriverPrecice.cc
  )

  if(PRECICE_LIBRARY)
    # preCICE built by the cfsdeps superbuild (cfsdeps/precice/External_PRECICE.cmake).
    # Consumed via link variables (the cfsdeps convention) - find_package cannot be
    # used here because the cfsdeps stage configures before the dependency is built.
    # The static libprecice pulls in boost / libxml2, which cfs links anyway.
    message(STATUS "preCICE from cfsdeps: ${PRECICE_LIBRARY}")
    # preciceadapter #includes precice headers at compile time, so it must be
    # built after the preCICE ExternalProject (which itself pulls in boost/eigen/libxml2).
    add_dependencies(preciceadapter precice)
    target_include_directories(preciceadapter PUBLIC
      ${PRECICE_INCLUDE_DIR} ${PRECICE_INCLUDE_DIR}/precice)
    target_link_libraries(preciceadapter PUBLIC
      ${PRECICE_LIBRARY} ${BOOST_LIBRARY} ${LIBXML2_LIBRARY})
  else()
    # Fallback: a preCICE installed on the system / provided by the developer
    # (e.g. a local dev build that does not use the cfsdeps superbuild).
    find_package(precice CONFIG REQUIRED)

    # Append the 'precice' subdirectory to the include paths
    # (the precice config misses the /precice path).
    get_target_property(precice_includes precice::precice INTERFACE_INCLUDE_DIRECTORIES)
    foreach(dir ${precice_includes})
      list(APPEND corrected_includes "${dir}/precice")
    endforeach()
    target_include_directories(preciceadapter PUBLIC ${corrected_includes})
    target_link_libraries(preciceadapter PUBLIC precice::precice)
  endif()
  #target_compile_definitions(preciceadapter PRIVATE USE_PRECICE)

  # Lightweight standalone preCICE "mock partner" used by the coupled testsuite
  # (Testsuite/TESTSUIT/preCICE). It replaces a heavy external solver (OpenFOAM)
  # in CI, so coupled tests can run cfs against deterministic analytic data.
  # Goes to ${EXECUTABLE_OUTPUT_PATH}=bin/ (and is packaged) like cfs/cfstool.
  add_executable(mock_fluid test/mock_fluid.cpp)
  if(PRECICE_LIBRARY)
    add_dependencies(mock_fluid precice)
    target_include_directories(mock_fluid PRIVATE
      ${PRECICE_INCLUDE_DIR} ${PRECICE_INCLUDE_DIR}/precice)
    target_link_libraries(mock_fluid PRIVATE
      ${PRECICE_LIBRARY} ${BOOST_LIBRARY} ${LIBXML2_LIBRARY})
  else()
    target_link_libraries(mock_fluid PRIVATE precice::precice)
  endif()
else()
  # Add DummyPreciceAdapter.cc to the preciceadapter target
  target_sources(preciceadapter PRIVATE
                 DummyPreciceAdapter.cc
                 IPreciceAdapter.cc
                 CoefFunctionGridElemDefaultPrecice.cc
                 CoefFunctionGridNodalDefaultPrecice.cc
                 TransientDriverPrecice.cc)
endif()

# Interface library for IPreciceAdapter
add_library(ipreciceadapter INTERFACE)
target_include_directories(ipreciceadapter INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})