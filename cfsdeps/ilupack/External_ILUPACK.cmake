# Ilupack is a closed source academic iterative solve by M. Bollhoefer, TU-Braunschweig
# by special agreement we have the source available which must not be redistributed

set(ilupack_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ilupack")
set(ilupack_source  "${ilupack_prefix}/src/ilupack")
set(ilupack_install  "${CMAKE_CURRENT_BINARY_DIR}")

# set(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -DPRINT_MEM -DPRINT_INLINE")
# set(CFSDEPS_ILUPACK_C_FLAGS "${CFSDEPS_C_FLAGS} -DPRINT_INLINE")
set(CFSDEPS_ILUPACK_C_FLAGS "${CFSDEPS_C_FLAGS}")

STRING(REPLACE ";" "^" ILUPACK_AMD_LIBRARY "${AMD_LIBRARY}")
STRING(REPLACE ";" "^" ILUPACK_BLAS_LIBRARY "${BLAS_LIBRARY}")
STRING(REPLACE ";" "^" ILUPACK_LAPACK_LIBRARY "${LAPACK_LIBRARY}")


#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
set(CMAKE_ARGS
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_INSTALL_PREFIX:PATH=${ilupack_install}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_ILUPACK_C_FLAGS}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DCMAKE_Fortran_COMPILER_ID:STRING=${CMAKE_Fortran_COMPILER_ID}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DAMD_LIB:STRING=${ILUPACK_AMD_LIBRARY}
  -DBLAS_LIB:FILEPATH=${ILUPACK_BLAS_LIBRARY}
  -DLAPACK_LIB:FILEPATH=${ILUPACK_LAPACK_LIBRARY}
  -DMETIS_LIB:FILEPATH=${METIS_LIBRARY}
  -DBUILD_SAMPLES:BOOL=OFF
)

if(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE})
endif()

if(CFS_DISTRO STREQUAL "MACOSX")
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT})
endif(CFS_DISTRO STREQUAL "MACOSX")

# decrypt source and patch it
set_from_env(CFS_KEY_ILUPACK)
set(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ilupack/ilupack-patch.cmake.in")
set(PFN "${ilupack_prefix}/ilupack-patch.cmake")
configure_file("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# The ilupack source is closed source! It must not be redistributed!!
#-------------------------------------------------------------------------------
set_from_env(CFS_DOWNLOAD_ILUPACK)
set(MIRRORS "${CFS_DS_SOURCES_DIR}/${CFS_DOWNLOAD_ILUPACK}/${ILUPACK_ZIP}")
set(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/ilupack/${ILUPACK_ZIP}")
set(MD5_SUM ${ILUPACK_MD5})

set(DLFN "${ilupack_prefix}/ilupack-download.cmake")
configure_file("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "ilupack" "${ILUPACK_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
set(TMP_DIR "${ilupack_prefix}")

set(ZIPFROMCACHE "${ilupack_prefix}/ilupackzipFromCache.cmake")
configure_file("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

set(ZIPTOCACHE "${ilupack_prefix}/ilupack-zipToCache.cmake")
configure_file("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of ILUPACK libraries.
#-----------------------------------------------------------------------------
set(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
set(ILUPACK_LIBRARY
  "${LD}/libDilupack.a;${LD}/libZilupack.a;${LD}/libblaslike.a;${LD}/libsparspak.a;${AMD_LIBRARY}" # AMD_LIBRARY already has the path
  CACHE FILEPATH "ILUPACK library.")
set(ILUPACK_DOUBLE_LIBRARY "${LD}/libDilupack.a")
set(ILUPACK_COMPLEX_LIBRARY "${LD}/libZilupack.a")
set(ILUPACK_BLASLIKE_LIBRARY "${LD}/libblaslike.a")
set(ILUPACK_SPARSPAK_LIBRARY "${LD}/libsparspak.a")
MARK_AS_ADVANCED(ILUPACK_LIBRARY)

#-------------------------------------------------------------------------------
# The ilupack external project
#-------------------------------------------------------------------------------
if("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(ilupack-double
    DEPENDS lapack metis suitesparse
    PREFIX "${ilupack_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "" )
  ExternalProject_Add(ilupack-complex
    PREFIX "${ilupack_prefix}"
    DOWNLOAD_COMMAND ""
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "" )

else()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project (double real)
  #-------------------------------------------------------------------------------
  
  # variables can be set in .cfs_platform_defaults.cmake, with cmake -D or as environment variables.
  if(NOT EXISTS LOCAL_FILE AND NOT DEFINED CFS_DOWNLOAD_ILUPACK)
    message(WARNING "no local ilupack sources exists and CFS_DOWNLOAD_ILUPACK not set.")
  endif()

  if(NOT DEFINED CFS_KEY_ILUPACK)
    message(WARNING "CFS_KEY_ILUPACK not set.")
  endif()  
  
  ExternalProject_Add(ilupack-double
    #DEPENDS ilupack_external_data lapack metis suitesparse
    DEPENDS lapack metis suitesparse
    PREFIX "${ilupack_prefix}"
    SOURCE_DIR "${ilupack_source}"
    URL ${LOCAL_FILE}
    URL_MD5 ${ILUPACK_MD5}
    # disable automatic extract of the zip file by cmake (does not work for encrypted files)
    DOWNLOAD_NO_EXTRACT 1
    # extract enceypted file and apply patches
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    LIST_SEPARATOR "^"
    CMAKE_ARGS ${CMAKE_ARGS} -DFLOAT_TYPE:STRING=DOUBLE_REAL
    BUILD_BYPRODUCTS ${ILUPACK_DOUBLE_LIBRARY} ${ILUPACK_BLASLIKE_LIBRARY} ${ILUPACK_SPARSPAK_LIBRARY}) 

  ExternalProject_Add(ilupack-complex
    DEPENDS ilupack-double
    PREFIX "${ilupack_prefix}"
    SOURCE_DIR "${ilupack_source}"
    DOWNLOAD_COMMAND ""
    LIST_SEPARATOR "^"
    CMAKE_ARGS ${CMAKE_ARGS} -DFLOAT_TYPE:STRING=DOUBLE_COMPLEX
    BUILD_BYPRODUCTS ${ILUPACK_COMPLEX_LIBRARY} )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(ilupack-double cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${ilupack_prefix} )

  ExternalProject_Add_Step(ilupack-complex cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${ilupack_prefix} )
  
  if("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache. complex seems sufficient
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(ilupack-complex cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR} )
  endif()
endif()

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
set(CFSDEPS ${CFSDEPS} ilupack-double ilupack-complex )

set(ILUPACK_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
mark_as_advanced(ILUPACK_INCLUDE_DIR)

