# Define the preciceadapter library
add_library(preciceadapter STATIC)

# Set the library's include directories to the current source directory
target_include_directories(preciceadapter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


if(USE_PRECICE)
  find_package(precice CONFIG REQUIRED)

  # Add PreciceAdapter.cc to the preciceadapter target
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

  # Retrieve the interface include directories from the imported target
  get_target_property(precice_includes precice::precice INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Using precice include directories: ${precice_includes}")

  # Append the 'precice' subdirectory to the existing include paths
  # I really don't know why precice config misses the /precice path...
  foreach(dir ${precice_includes})
    list(APPEND corrected_includes "${dir}/precice")
  endforeach()
  message(STATUS "Corrected precice include directories: ${corrected_includes}")
  # Add these include directories to the target
  target_include_directories(preciceadapter PUBLIC ${corrected_includes})
  target_link_libraries(preciceadapter PUBLIC precice::precice)
  #target_compile_definitions(preciceadapter PRIVATE USE_PRECICE)
else()
  # Add DummyPreciceAdapter.cc to the preciceadapter target
  target_sources(preciceadapter PRIVATE 
                 DummyPreciceAdapter.cc
                 IPreciceAdapter.cc
                 CoefFunctionGridElemDefaultPrecice.cc
                 TransientDriverPrecice.cc)
endif()

# Interface library for IPreciceAdapter
add_library(ipreciceadapter INTERFACE)
target_include_directories(ipreciceadapter INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})