#-------------------------------------------------------------------------------
# Set prefix path and path to hdf5 sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(hdf5_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/hdf5")
set(hdf5_install  "${CMAKE_CURRENT_BINARY_DIR}")
set(hdf5_source  "${hdf5_prefix}/src/hdf5")

#SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/hdf5/hdf5-patch.cmake.in")
#SET(PFN "${hdf5_prefix}/hdf5-patch.cmake")
#CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${hdf5_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
#  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${hdf5_install}/bin/${CFS_ARCH_STR}
#  -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${hdf5_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}
#  -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${hdf5_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}
  -DHDF5_INSTALL_BIN_DIR:PATH=bin/${CFS_ARCH_STR}
  -DHDF5_INSTALL_LIB_DIR:PATH=${LIB_SUFFIX}/${CFS_ARCH_STR}
  # We do not want to see warning messages from external projects
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
)

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

#-------------------------------------------------------------------------------
# The hdf5-static external project
# We do not need it at the moment since we use the HDF5 inside ParaView for
# our plugin.
#-------------------------------------------------------------------------------
#ExternalProject_Add(hdf5-shared
#  DEPENDS zlib
#  PREFIX ${hdf5_prefix}
#  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/hdf5
#  SOURCE_DIR ${hdf5_source}
#  URL ${HDF5_URL}/${HDF5_GZ}
#  URL_MD5 ${HDF5_MD5}
#  CMAKE_ARGS
#     ${CMAKE_ARGS}
#    -DBUILD_SHARED_LIBS:BOOL=ON
#    -DHDF5_BUILD_CPP_LIB:BOOL=ON
#    -DHDF5_BUILD_HL_LIB:BOOL=ON
#    -DHDF5_BUILD_FORTRAN:BOOL=OFF
#    -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
#    -DZLIB_INCLUDE_DIR:PATH=${hdf5_install}/include
#    -DZLIB_LIBRARY:FILEPATH=${hdf5_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_SHARED_LIBRARY_PREFIX}z${CMAKE_SHARED_LIBRARY_SUFFIX}
#    )

IF(MINGW)
  IF(CMAKE_CROSSCOMPILING)
    LIST(APPEND CMAKE_ARGS
      -DH5_HAVE_GETHOSTNAME:INTERNAL=0
      )
    
    IF(CFS_TARGET_ARCH STREQUAL "X86_64")
      LIST(APPEND CMAKE_ARGS
	-C ${CMAKE_CURRENT_SOURCE_DIR}/cfsdeps/hdf5/TryRunResults_MINGW_X86_64_CENTOS6.cmake
	)
    ENDIF()
  ENDIF()

ENDIF(MINGW)

IF(APPLE)
  IF(CMAKE_CROSSCOMPILING)
    IF(CFS_TARGET_ARCH STREQUAL "X86_64")
      LIST(APPEND CMAKE_ARGS
	-C ${CMAKE_CURRENT_SOURCE_DIR}/cfsdeps/hdf5/TryRunResults_MACOSX_X86_64_CENTOS6.cmake
	)
    ELSE()
      LIST(APPEND CMAKE_ARGS
	-C ${CMAKE_CURRENT_SOURCE_DIR}/cfsdeps/hdf5/TryRunResults_MACOSX_I386_CENTOS6.cmake
	)
    ENDIF()
  ENDIF()
ENDIF(APPLE)


LIST(APPEND CMAKE_ARGS
  -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
  )

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/hdf5/hdf5-patch.cmake.in")
SET(PFN "${hdf5_prefix}/hdf5-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since lse17 may not be
# accessible from behind firewalls.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.8.tar.gz"
  "ftp://ftp.ca.freebsd.org/pub/packages/graphics/netcdf/netcdf-4/hdf5-1.8.8.tar.gz"
#  "http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.8/src/hdf5-1.8.8.tar.gz"
)

#-------------------------------------------------------------------------------
# Try to download sources into CFSDEPS cache directory.
#-------------------------------------------------------------------------------
DOWNLOAD_CFSDEPS(
  "${CFS_DEPS_CACHE_DIR}/sources/hdf5/${HDF5_GZ}"
  ${HDF5_MD5}
  "${MIRRORS}"
)


#-------------------------------------------------------------------------------
# The hdf5-static external project
#-------------------------------------------------------------------------------
ExternalProject_Add(hdf5-static
  DEPENDS zlib
  PREFIX ${hdf5_prefix}
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/hdf5
  SOURCE_DIR ${hdf5_source}
  URL ${HDF5_URL}/${HDF5_GZ}
  URL_MD5 ${HDF5_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  LIST_SEPARATOR ^
  CMAKE_ARGS
    ${CMAKE_ARGS}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DHDF5_BUILD_CPP_LIB:BOOL=ON
    -DHDF5_BUILD_HL_LIB:BOOL=ON
    -DHDF5_BUILD_FORTRAN:BOOL=OFF
    -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
    -DZLIB_INCLUDE_DIR:PATH=${hdf5_install}/include
    -DHDF5_BUILD_TOOLS:BOOL=ON
    # On Mac OS X we can get problems with the system strdup function.
    -DH5_HAVE_STRDUP:BOOL=OFF
    )

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  hdf5-static
)

#-------------------------------------------------------------------------------
# Determine paths of HDF5 libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  IF(MINGW)
    SET(HDF5_LIBRARY_DEBUG
      "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
    SET(HDF5_LIBRARY_RELEASE
      "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5${CMAKE_STATIC_LIBRARY_SUFFIX}")

    SET(HDF5_CPP_LIBRARY_DEBUG
      "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_cpp_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
    SET(HDF5_CPP_LIBRARY_RELEASE
      "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_cpp${CMAKE_STATIC_LIBRARY_SUFFIX}")

    SET(HDF5_LT_LIBRARY_DEBUG
      "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
    SET(HDF5_LT_LIBRARY_RELEASE
      "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl${CMAKE_STATIC_LIBRARY_SUFFIX}")

    SET(HDF5_LT_CPP_LIBRARY_DEBUG
      "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_cpp_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
    SET(HDF5_LT_CPP_LIBRARY_RELEASE
      "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_cpp${CMAKE_STATIC_LIBRARY_SUFFIX}")
  ELSE()
    SET(HDF5_LIBRARY_RELEASE ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/hdf5dll.lib)
    SET(HDF5_CPP_LIBRARY_RELEASE ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/hdf5_cppdll.lib)
    SET(HDF5_LIBRARY_DEBUG ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/hdf5ddll.lib)
    SET(HDF5_CPP_LIBRARY_DEBUG ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/hdf5_cppddll.lib)
  ENDIF()
ELSE(WIN32)
  SET(HDF5_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5${CMAKE_STATIC_LIBRARY_SUFFIX}")

  SET(HDF5_CPP_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_cpp_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_CPP_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_cpp${CMAKE_STATIC_LIBRARY_SUFFIX}")

  SET(HDF5_LT_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_LT_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl${CMAKE_STATIC_LIBRARY_SUFFIX}")

  SET(HDF5_LT_CPP_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_cpp_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_LT_CPP_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_cpp${CMAKE_STATIC_LIBRARY_SUFFIX}")

ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Set HDF5_LIBRARY* according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(HDF5_LIBRARY "${HDF5_LIBRARY_DEBUG}")
  SET(HDF5_CPP_LIBRARY "${HDF5_CPP_LIBRARY_DEBUG}")

  SET(HDF5_LT_LIBRARY "${HDF5_LT_LIBRARY_DEBUG}")
  SET(HDF5_LT_CPP_LIBRARY "${HDF5_LT_CPP_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(HDF5_LIBRARY "${HDF5_LIBRARY_RELEASE}")
  SET(HDF5_CPP_LIBRARY "${HDF5_CPP_LIBRARY_RELEASE}")

  SET(HDF5_LT_LIBRARY "${HDF5_LT_LIBRARY_RELEASE}")
  SET(HDF5_LT_CPP_LIBRARY "${HDF5_LT_CPP_LIBRARY_RELEASE}")
ENDIF(DEBUG)

SET(HDF5_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "hdf5 include directory")

INCLUDE_DIRECTORIES("${CFS_BINARY_DIR}/include" "${CFS_BINARY_DIR}/include/cpp")

MARK_AS_ADVANCED(HDF5_INCLUDE_DIR)
