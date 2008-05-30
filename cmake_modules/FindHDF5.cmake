SET(HDF5_FOUND 0)

#-------------------------------------------------------------------------------
# Find Zlib.
#-------------------------------------------------------------------------------
FIND_PACKAGE(ZLIB)

MARK_AS_ADVANCED(ZLIB_LIBRARY)
MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# Determine paths of HDF5 libraries.
#-------------------------------------------------------------------------------
FIND_LIBRARY(HDF5_LIBRARY
  NAMES hdf5
  PATHS ${CFSDEPS_LIBRARY_DIR}
  )

FIND_LIBRARY(HDF5_CPP_LIBRARY
  NAMES hdf5_cpp
  PATHS ${CFSDEPS_LIBRARY_DIR}
  )

#-------------------------------------------------------------------------------
# Mark paths of HDF5 libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(HDF5_LIBRARY)
MARK_AS_ADVANCED(HDF5_CPP_LIBRARY)

FIND_PATH(HDF5_INCLUDE_DIR hdf5.h 
  ${CFSDEPS_INCLUDE_DIR} 
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH)

MARK_AS_ADVANCED(HDF5_INCLUDE_DIR)


IF(HDF5_LIBRARY AND HDF5_CPP_LIBRARY AND HDF5_INCLUDE_DIR)
  SET(HDF5_FOUND 1)
ENDIF(HDF5_LIBRARY AND HDF5_CPP_LIBRARY AND HDF5_INCLUDE_DIR)


IF(HDF5_FOUND)

  #-----------------------------------------------------------------------------
  # Determine version of HDF5, zlib and szlib by
  # preprocessing their version header
  #-----------------------------------------------------------------------------
  SET(HDF5_TEST_FILE "${CFS_BINARY_DIR}/include/get_hdf5_infos.hh")

  FILE(WRITE ${HDF5_TEST_FILE}
    "#include <${HDF5_INCLUDE_DIR}/H5public.h>\n\n"
    "#include <${ZLIB_INCLUDE_DIR}/zlib.h>\n\n"
    "<CFS_HDF5_VERSION>H5_VERS_INFO</CFS_HDF5_VERSION>\n"
    "<CFS_ZLIB_VERSION>ZLIB_VERSION</CFS_ZLIB_VERSION>\n")

  SET(HDF5_VERSION_CMD "${CMAKE_CXX_COMPILER} -E ${HDF5_TEST_FILE}")

  EXEC_PROGRAM(${HDF5_VERSION_CMD}
    ARGS
    OUTPUT_VARIABLE HDF5_VERSION_STR
    RETURN_VALUE RETVAL)

#  MESSAGE("HDF5_VERSION_STR ${HDF5_VERSION_STR}")
  FILE(REMOVE ${HDF5_TEST_FILE})

  STRING(REGEX MATCH
    "<CFS_HDF5_VERSION>.*</CFS_HDF5_VERSION>"
    CFS_HDF5_VERSION
    "${HDF5_VERSION_STR}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
        CFS_HDF5_VERSION "${CFS_HDF5_VERSION}")

  STRING(REGEX MATCH
    "<CFS_ZLIB_VERSION>.*</CFS_ZLIB_VERSION>"
    CFS_ZLIB_VERSION
    "${HDF5_VERSION_STR}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
        CFS_ZLIB_VERSION "${CFS_ZLIB_VERSION}")

#  MESSAGE("CFS_HDF5_VERSION ${CFS_HDF5_VERSION}")
#  MESSAGE("CFS_ZLIB_VERSION ${CFS_ZLIB_VERSION}")

  STRING(REGEX REPLACE "hdf5\\.so" "hdf5.a;-lpthread;-ldl"
    HDF5_LIBRARY "${HDF5_LIBRARY}")

  SET(HDF5_LIBRARY "${HDF5_LIBRARY}" CACHE PATH "Path to HDF5 libraries." FORCE)

ENDIF(HDF5_FOUND)

IF(NOT HDF5_FOUND)
  MESSAGE("Warning: HDF5 could not be found! Please specify proper paths.")
ENDIF(NOT HDF5_FOUND)
