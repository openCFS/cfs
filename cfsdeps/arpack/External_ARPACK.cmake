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

#-------------------------------------------------------------------------------
# The ARPACK external project
#-------------------------------------------------------------------------------
ExternalProject_Add(arpack
  PREFIX "${ARPACK_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/arpack
  URL ${ARPACK_URL}/${ARPACK_GZ}
  URL_MD5 ${ARPACK_MD5}
  CMAKE_ARGS
    ${CMAKE_ARGS}
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/arpack/arpack-patch.cmake.in")
SET(PFN "${ARPACK_prefix}/arpack-patch.cmake")
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
ExternalProject_Add_Step(arpack custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${ARPACK_source}
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
