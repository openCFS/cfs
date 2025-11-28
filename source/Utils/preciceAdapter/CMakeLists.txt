# Define the preciceadapter library
add_library(preciceadapter STATIC
  IPreciceAdapter.cc
  DummyPreciceAdapter.cc
  CoefFunctionGridElemDefaultPrecice.cc
  CoefFunctionGridNodalDefaultPrecice.cc
  TransientDriverPrecice.cc
)

# Set the library's include directories to the current source directory
target_include_directories(preciceadapter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(USE_PRECICE)
  target_sources(preciceadapter PRIVATE 
    MinimalXmlParser.cc
    PreciceAdapter.cc
    NodeResult.cc
    ElementResult.cc
  )

  if(NOT TARGET precice::precice)
    message(FATAL_ERROR "USE_PRECICE is ON but precice::precice target is missing. Check External_PRECICE.cmake.")
  endif()

  # Retrieve the interface include directories from the imported target
  get_target_property(precice_includes precice::precice INTERFACE_INCLUDE_DIRECTORIES)
  if(precice_includes)
    message(STATUS "Using precice include directories: ${precice_includes}")

    # Append the 'precice' subdirectory to the existing include paths
    # I really don't know why precice config misses the /precice path...
    set(corrected_includes "")
    foreach(dir ${precice_includes})
      list(APPEND corrected_includes "${dir}/precice")
    endforeach()
    message(STATUS "Corrected precice include directories: ${corrected_includes}")
    # Add these include directories to the target
    target_include_directories(preciceadapter PUBLIC ${corrected_includes})
  endif()
  target_link_libraries(preciceadapter PUBLIC precice::precice)
endif()

# Interface library for IPreciceAdapter
add_library(ipreciceadapter INTERFACE)
target_include_directories(ipreciceadapter INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
