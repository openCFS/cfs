#-------------------------------------------------------------------------------
# CFD General Notation System (CGNS)
# Needed for ADF routines by STARCCM+ reader. Also provides adfview.
#                                           
# Project Homepage                          
# http://www.cgns.org                       
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to CGNS sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(cgns_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cgns")
set(cgns_install  "${CMAKE_CURRENT_BINARY_DIR}")
set(cgns_source  "${cgns_prefix}/src/cgns")

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${cgns_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DEXECUTABLE_OUTPUT_PATH:PATH=${cgns_install}/bin/${CFS_ARCH_STR}
  -DLIBRARY_OUTPUT_PATH:PATH=${cgns_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DBUILD_CGNSTOOLS:BOOL=ON
  -DENABLE_HDF5:BOOL=ON
  -DHDF5_INCLUDE_PATH:PATH=${cgns_install}/include
  -DHDF5_LIBRARY:FILEPATH=${HDF5_LIBRARY},-lm
  -DHDF5_NEED_ZLIB:BOOL=ON
  -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
  -DCGNS_BUILD_SHARED:BOOL=OFF
  -DCGNS_USE_SHARED:BOOL=OFF
  -DBUILD_SHARED_LIBS:BOOL=OFF
  # We do not want to see warning messages from external projects
  -DCMAKE_C_FLAGS:STRING=${CFLAGS}
  -DCMAKE_CXX_FLAGS:STRING=${CFLAGS}
)

IF(CFS_DISTRO STREQUAL "MACOSX")
  # Explicitly set build architectures and  system SDK root dir to match
  # those given in CMake.
  
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    "-DCMAKE_EXE_LINKER_FLAGS:STRING=-L/usr/X11/lib -lX11 -lXmu"
  )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgns/cgns-patch.cmake.in")
SET(PFN "${cgns_prefix}/cgns-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://heanet.dl.sourceforge.net/project/cgns/cgnslib_3.1/cgnslib_3.1-2.tar.gz"
  "http://mirror.transact.net.au/pub/sourceforge/c/project/cg/cgns/cgnslib_3.1/cgnslib_3.1.3-2.tar.gz"
  "${CGNS_URL}/${CGNS_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/cgns/${CGNS_GZ}")
SET(MD5_SUM ${CGNS_MD5})

SET(DLFN "${cgns_prefix}/cgns-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
  )

#-------------------------------------------------------------------------------
# The CGNS-static external project
#-------------------------------------------------------------------------------
ExternalProject_Add(cgns-static
  DEPENDS hdf5-static zlib
  PREFIX ${cgns_prefix}
  SOURCE_DIR ${cgns_source}
  URL ${LOCAL_FILE}
  URL_MD5 ${CGNS_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  LIST_SEPARATOR ,
  CMAKE_ARGS
     ${CMAKE_ARGS}
    )

#-------------------------------------------------------------------------------
# Add custom download step to be able to download from a list of mirrors
# instead of just a single URL.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(cgns-static cfsdeps_download
   COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
   DEPENDERS download
   DEPENDS "${DLFN}"
   WORKING_DIRECTORY ${cgns_prefix}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  cgns-static
)

#-------------------------------------------------------------------------------
# Determine paths of CGNS libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(CGNS_LIBRARY_RELEASE ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/CGNSdll.lib)
  SET(CGNS_CPP_LIBRARY_RELEASE ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/cgns_cppdll.lib)
  SET(CGNS_LIBRARY_DEBUG ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/CGNSddll.lib)
  SET(CGNS_CPP_LIBRARY_DEBUG ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/cgns_cppddll.lib)
ELSE(WIN32)
  SET(CGNS_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(CGNS_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}CGNS${CMAKE_STATIC_LIBRARY_SUFFIX}")

  SET(CGNS_CPP_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns_cpp_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(CGNS_CPP_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns_cpp${CMAKE_STATIC_LIBRARY_SUFFIX}")

  SET(CGNS_LT_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns_hl_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(CGNS_LT_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns_hl${CMAKE_STATIC_LIBRARY_SUFFIX}")

  SET(CGNS_LT_CPP_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns_hl_cpp_debug${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(CGNS_LT_CPP_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns_hl_cpp${CMAKE_STATIC_LIBRARY_SUFFIX}")

ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Set CGNS_LIBRARY* according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(CGNS_LIBRARY "${CGNS_LIBRARY_DEBUG}")
  SET(CGNS_CPP_LIBRARY "${CGNS_CPP_LIBRARY_DEBUG}")

  SET(CGNS_LT_LIBRARY "${CGNS_LT_LIBRARY_DEBUG}")
  SET(CGNS_LT_CPP_LIBRARY "${CGNS_LT_CPP_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(CGNS_LIBRARY "${CGNS_LIBRARY_RELEASE}")
  SET(CGNS_CPP_LIBRARY "${CGNS_CPP_LIBRARY_RELEASE}")

  SET(CGNS_LT_LIBRARY "${CGNS_LT_LIBRARY_RELEASE}")
  SET(CGNS_LT_CPP_LIBRARY "${CGNS_LT_CPP_LIBRARY_RELEASE}")
ENDIF(DEBUG)

SET(CGNS_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "CGNS include directory")

INCLUDE_DIRECTORIES("${CFS_BINARY_DIR}/include" "${CFS_BINARY_DIR}/include/cpp")

MARK_AS_ADVANCED(CGNS_INCLUDE_DIR)
