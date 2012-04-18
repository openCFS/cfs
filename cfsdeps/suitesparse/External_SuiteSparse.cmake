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
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} library
  INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
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

SET(CHOLMOD_INCLUDE_DIR "${CFSDEPS_INCLUDE_DIR}")

#-------------------------------------------------------------------------------
# Determine paths of CholMod libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(CHOLMOD_LIBRARY
  "${LD}/libcholmod.a;${LD}/libamd.a;${LD}/libcolamd.a;${LD}/libcamd.a;${LD}/libccolamd.a"
  CACHE FILEPATH "CholMod library.")

MARK_AS_ADVANCED(CHOLMOD_LIBRARY)

