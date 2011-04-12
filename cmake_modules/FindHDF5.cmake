SET(HDF5_FOUND 0)

#-------------------------------------------------------------------------------
# Look for HDF5 header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("HDF5"
  "${CFS_BINARY_DIR}/include/hdf5.h"
  "${CFS_DEPS_ROOT}/hdf5/build_hdf5.pl"
  "build_hdf5.log")

#-------------------------------------------------------------------------------
# Determine paths of HDF5 libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(HDF5_LIBRARY_RELEASE ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/hdf5dll.lib)
  SET(HDF5_CPP_LIBRARY_RELEASE ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/hdf5_cppdll.lib)
  SET(HDF5_LIBRARY_DEBUG ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/hdf5ddll.lib)
  SET(HDF5_CPP_LIBRARY_DEBUG ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/hdf5_cppddll.lib)
ELSE(WIN32)
  SET(HDF5_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libhdf5.a")

  SET(HDF5_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libhdf5.a")

  SET(HDF5_CPP_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libhdf5_cpp.a")
    
  SET(HDF5_CPP_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libhdf5_cpp.a")
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Set HDF5_LIBRARY* according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(HDF5_LIBRARY "${HDF5_LIBRARY_DEBUG}")
  SET(HDF5_CPP_LIBRARY "${HDF5_CPP_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(HDF5_LIBRARY "${HDF5_LIBRARY_RELEASE}")
  SET(HDF5_CPP_LIBRARY "${HDF5_CPP_LIBRARY_RELEASE}")
ENDIF(DEBUG)

SET(HDF5_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

# Safety check for existance of libs and headers
IF(EXISTS "${HDF5_LIBRARY}"
    AND EXISTS "${HDF5_CPP_LIBRARY}"
    AND EXISTS "${HDF5_INCLUDE_DIR}/hdf5.h")
  
  SET(HDF5_FOUND 1)
  INCLUDE_DIRECTORIES("${CFS_BINARY_DIR}/include" "${CFS_BINARY_DIR}/include/cpp")
  
ELSE(EXISTS "${HDF5_LIBRARY}"
    AND EXISTS "${HDF5_CPP_LIBRARY}"
    AND EXISTS "${HDF5_INCLUDE_DIR}/hdf5.h")
  
  MESSAGE(FATAL_ERROR "Could not find ${HDF5_LIBRARY}, ${HDF5_CPP_LIBRARY} or"
    " ${CFS_BINARY_DIR}/include/hdf5.h")
  
ENDIF(EXISTS "${HDF5_LIBRARY}"
  AND EXISTS "${HDF5_CPP_LIBRARY}"
  AND EXISTS "${HDF5_INCLUDE_DIR}/hdf5.h")


IF(NOT WIN32)
SET(HDF5_LIBRARY "${HDF5_LIBRARY};${ZLIB_LIBRARY};-lpthread")
ENDIF(NOT WIN32)

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
  SET(HDF5_CPP_LIBRARY "${HDF5_CPP_LIBRARY}" CACHE PATH "Path to HDF5 C++ libraries." FORCE)

  MARK_AS_ADVANCED(HDF5_CPP_LIBRARY)
  MARK_AS_ADVANCED(HDF5_LIBRARY)
ENDIF(HDF5_FOUND)

# MESSAGE("CFS_HDF5_VERSION ${CFS_HDF5_VERSION}")

IF(WIN32)
  IF(DEBUG)
    SET(HDF5DLL "${CFS_BINARY_DIR}/bin/${CFS_ARCH_STR}/hdf5ddll.dll")
    SET(HDF5_CPPDLL "${CFS_BINARY_DIR}/bin/${CFS_ARCH_STR}/hdf5_cppddll.dll")
  ELSE(DEBUG)
    SET(HDF5DLL "${CFS_BINARY_DIR}/bin/${CFS_ARCH_STR}/hdf5dll.dll")
    SET(HDF5_CPPDLL "${CFS_BINARY_DIR}/bin/${CFS_ARCH_STR}/hdf5_cppdll.dll")
  ENDIF(DEBUG)
ENDIF(WIN32)

IF(NOT HDF5_FOUND)
  MESSAGE("Warning: HDF5 could not be found! Please specify proper paths.")
ENDIF(NOT HDF5_FOUND)
