#-------------------------------------------------------------------------------
# SuiteSparse:
#
#      SuiteSparse is a single archive that contains all packages that I have
#      authored or co-authored that are available at this site. This gives you
#      a simple way of getting and installing all of my software packages.
#      Currently, this includes:
#          o AMD: symmetric approximate minimum degree
#          o BTF: permutation to block triangular form
#          o CAMD: symmetric approximate minimum degree
#          o CCOLAMD: constrained column approximate minimum degree
#          o COLAMD: column approximate minimum degree
#          o CHOLMOD: sparse supernodal Cholesky factorization and update/downdate
#          o CSparse: a concise sparse matrix package
#          o CXSparse: an extended version of CSparse
#          o KLU: sparse LU factorization, for circuit simulation
#          o LDL: a simple LDL^T factorization
#          o UMFPACK: sparse multifrontal LU factorization
#          o RBio: MATLAB toolbox for reading/writing sparse matrices
#          o UFconfig: common configuration for all but CSparse
#          o LINFACTOR: solve Ax=b using LU or CHOL
#          o MESHND: 2D and 3D mesh generation and nested dissection
#          o SSMULT: sparse matrix times sparse matrix
#          o SuiteSparseQR: multifrontal sparse QR 
#
# Project Homepage
#
# http://www.cise.ufl.edu/research/sparse/SuiteSparse/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to suitesparse sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(suitesparse_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/suitesparse")
set(suitesparse_source  "${suitesparse_prefix}/src/suitesparse")
set(suitesparse_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${suitesparse_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCMAKE_LINKER:FILEPATH=${CMAKE_LINKER}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DCFS_FORTRAN_COMPILER_NAME=${CFS_FORTRAN_COMPILER_NAME}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCFS_INCLUDE_DIR:PATH=${CFS_BINARY_DIR}/include
  )


IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/suitesparse/suitesparse-patch.cmake.in")
SET(PFN "${suitesparse_prefix}/suitesparse-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://www.cise.ufl.edu/research/sparse/SuiteSparse/${SUITESPARSE_GZ}"
  "${SUITESPARSE_URL}/${SUITESPARSE_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/suitesparse/${SUITESPARSE_GZ}")
SET(MD5_SUM ${SUITESPARSE_MD5})

SET(DLFN "${suitesparse_prefix}/suitesparse-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
) 

IF(WIN32)
  SET(PRECOMPILED_PCKG_NAME "suitesparse_${SUITESPARSE_VER}_${CFS_ARCH_STR}_${TOOLSET_ID}_${CMAKE_BUILD_TYPE}.zip")
ELSE(WIN32)
  SET(PRECOMPILED_PCKG_NAME "suitesparse_${SUITESPARSE_VER}_${CFS_ARCH_STR}_${FC_ID}_${CMAKE_BUILD_TYPE}.zip")
ENDIF(WIN32)
SET(PRECOMPILED_PCKG_FILE "${CFS_DEPS_CACHE_DIR}/precompiled/CFSDEPS/${PRECOMPILED_PCKG_NAME}")
  
SET(PREFIX_DIR "${suitesparse_prefix}")

SET(ZIPFROMCACHE "${suitesparse_prefix}/suitesparse-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${suitesparse_prefix}/suitesparse-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The suitesparse external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(suitesparse
    PREFIX "${suitesparse_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(suitesparse
    DEPENDS metis
    PREFIX "${suitesparse_prefix}"
    SOURCE_DIR "${suitesparse_source}"
    URL ${LOCAL_FILE}
    URL_MD5 ${SUITESPARSE_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=OFF
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(suitesparse cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${suitesparse_prefix}
  )
  
  IF("${CFS_DEPS_TOCACHE}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(suitesparse cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  suitesparse
)

SET(CHOLMOD_INCLUDE_DIR "${CFSDEPS_INCLUDE_DIR}")

#-------------------------------------------------------------------------------
# Determine paths of CholMod libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")

SET(AMD_LIBS
  amd_dint
  amd_dlong
  amd)

SET(AMD_LIBRARY "")

foreach(lib IN LISTS AMD_LIBS)
  LIST(APPEND AMD_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}")
endforeach()

SET(CHOLMOD_LIBS
  cholmod_dlong
  cholmod_dint
  colamd_dint
  colamd_dlong
  colamd
  camd_dint
  camd_dlong
  camd
  ccolamd_dint
  ccolamd_dlong
  ccolamd
  SuiteSparse_config
  )

SET(CHOLMOD_LIBRARY "")

foreach(lib IN LISTS CHOLMOD_LIBS)
  LIST(APPEND CHOLMOD_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}")
endforeach()

IF(NOT CFS_DISTRO STREQUAL "MACOSX")
  IF(NOT MINGW)
    LIST(APPEND CHOLMOD_LIBRARY
      "-lrt"
      )
  ENDIF()
ENDIF()

LIST(APPEND CHOLMOD_LIBRARY
  ${AMD_LIBRARY}
  ${METIS_LIBRARY}
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY})

SET(UMFPACK_LIBS
  umfpack_dlong
  umfpack_dint
  umfpack_zlong
  umfpack_zint
  umfpack
  )

SET(UMFPACK_LIBRARY "")

foreach(lib IN LISTS UMFPACK_LIBS)
  LIST(APPEND UMFPACK_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}")
endforeach()

LIST(APPEND UMFPACK_LIBRARY
  ${CHOLMOD_LIBRARY})
