IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} LESS 2.8 AND
    NOT EXISTS "${CFS_BINARY_DIR}/cmake")
#    MESSAGE ("It is less than 2.8")

  #-----------------------------------------------------------------------------
  # HDF5 and ParaView 3.10 require CMake 2.8
  #-----------------------------------------------------------------------------
  BUILD_EXTLIB("CMake 2.8.4"
    "${CFS_BINARY_DIR}/cmake"
    "${CFS_DEPS_ROOT}/cmake/build_cmake.pl"
    "build_cmake.log")
ENDIF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} LESS 2.8 AND
        NOT EXISTS "${CFS_BINARY_DIR}/cmake")
