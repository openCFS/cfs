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
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCFS_INCLUDE_DIR:PATH=${CFS_BINARY_DIR}/include
  )


IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_C_FLAGS:FILEPATH=${CFLAGS}
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
# The suitesparse external project
#-------------------------------------------------------------------------------
ExternalProject_Add(suitesparse
  DEPENDS metis
  PREFIX "${suitesparse_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/suitesparse
  SOURCE_DIR "${suitesparse_source}"
  URL ${SUITESPARSE_URL}/${SUITESPARSE_GZ}
  URL_MD5 ${SUITESPARSE_MD5}
  CMAKE_ARGS
    ${CMAKE_ARGS}
    -DBUILD_SHARED_LIBS:BOOL=OFF
)


#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/suitesparse/suitesparse-patch.cmake.in")
SET(PFN "${suitesparse_prefix}/suitesparse-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# We do not use the PATCH_COMMAND  of ExternalProject_Add since we do not only
# want to apply the patch script  during configuration time but also if it has
# changed.  Therefore,   we  need  a   dependency  on  the   configured  patch
# script. This can be achieved by  adding an additional build step between the
# download and configure steps.
#
# NOTE: The  patch script should  be designed  in such a  way, that it  can be
# applied to  an already patched  source tree. This  is due to the  fact, that
# ExternalProject_Add only extracts the source if the MD5 sum has has changed.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(suitesparse custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${suitesparse_source}
)

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
  )

SET(CHOLMOD_LIBRARY "")

foreach(lib IN LISTS CHOLMOD_LIBS)
  LIST(APPEND CHOLMOD_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}")
endforeach()

LIST(APPEND CHOLMOD_LIBRARY
  ${AMD_LIBRARY}
  ${METIS_LIBRARY})

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
