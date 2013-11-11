#-------------------------------------------------------------------------------
# ARnoldi PACKage -  ARPACK is a collection of  Fortran77 subroutines designed
# to solve large scale eigenvalue problems.
#
# Project Homepage
# http://www.caam.rice.edu/software/ARPACK/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to ARPACK sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(ARPACK_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/arpack")
set(ARPACK_source  "${ARPACK_prefix}/src/arpack")
set(ARPACK_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${ARPACK_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=-w
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
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
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/arpack/arpack-patch.cmake.in")
SET(PFN "${ARPACK_prefix}/arpack-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://ftp1.rrzn.uni-hannover.de/pub/mirror/bsd/FreeBSD/ports/distfiles/arpack-ng_3.1.1.tar.gz"
  "ftp://ftp.uwsg.indiana.edu/pub/FreeBSD/ports/distfiles/arpack-ng_3.1.1.tar.gz"
  "${ARPACK_URL}/${ARPACK_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/arpack/${ARPACK_GZ}")
SET(MD5_SUM ${ARPACK_MD5})

SET(DLFN "${ARPACK_prefix}/arpack-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
  ) 

#-------------------------------------------------------------------------------
# The ARPACK external project
#-------------------------------------------------------------------------------
ExternalProject_Add(arpack
  PREFIX "${ARPACK_prefix}"
  URL ${LOCAL_FILE}
  URL_MD5 ${ARPACK_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CMAKE_ARGS
    ${CMAKE_ARGS}
)

#-------------------------------------------------------------------------------
# Add custom download step to be able to download from a list of mirrors
# instead of just a single URL.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(arpack cfsdeps_download
   COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
   DEPENDERS download
   DEPENDS "${DLFN}"
   WORKING_DIRECTORY ${ARPACK_prefix}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  arpack
)

#-------------------------------------------------------------------------------
# Determine paths of ARPACK libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(ARPACK_LIBRARY
  "${LD}/libarpack.a"
  CACHE FILEPATH "ARPACK library.")

# MARK_AS_ADVANCED(ARPACK_INCLUDE_DIR)
MARK_AS_ADVANCED(ARPACK_LIBRARY)
