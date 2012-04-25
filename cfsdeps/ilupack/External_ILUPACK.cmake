#-------------------------------------------------------------------------------
# Set paths to ilupack sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(ilupack_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ilupack")
set(ilupack_source  "${ilupack_prefix}/src/ilupack")
set(ilupack_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_INSTALL_PREFIX:PATH=${ilupack_install}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=-w
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=-w
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DCFS_DISTRO:STRING=${CFS_DISTRO}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
)

#-------------------------------------------------------------------------------
# The ilupack external project (double real)
#-------------------------------------------------------------------------------
ExternalProject_Add(ilupack-double
  DEPENDS lapack metis suitesparse
  PREFIX "${ilupack_prefix}"
  SOURCE_DIR "${ilupack_source}"
  URL ${ILUPACK_PATH}/${ILUPACK_GZ}
  URL_MD5 ${ILUPACK_MD5}
  CMAKE_ARGS
    ${CMAKE_ARGS}
    -DFLOAT_TYPE:STRING=DOUBLE_REAL
    -DBUILD_SAMPLES:BOOL=
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ilupack/ilupack-patch.cmake.in")
SET(PFN "${ilupack_prefix}/ilupack-patch.cmake")
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
ExternalProject_Add_Step(ilupack-double custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${ilupack_source}
)

#-------------------------------------------------------------------------------
# The ilupack external project (double complex)
#-------------------------------------------------------------------------------
ExternalProject_Add(ilupack-complex
  DEPENDS ilupack-double
  PREFIX "${ilupack_prefix}"
  SOURCE_DIR "${ilupack_source}"
  DOWNLOAD_COMMAND ""
  CMAKE_ARGS
    ${CMAKE_ARGS}
    -DFLOAT_TYPE:STRING=DOUBLE_COMPLEX
    -DBUILD_SAMPLES:BOOL=
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  ilupack-double
  ilupack-complex
)

SET(ILUPACK_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of ILUPACK libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(ILUPACK_LIBRARY
  "${LD}/libDilupack.a;${LD}/libZilupack.a;${LD}/libblaslike.a;${LD}/libsparspak.a;${LD}/libamd.a"
  CACHE FILEPATH "ILUPACK library.")

MARK_AS_ADVANCED(ILUPACK_LIBRARY)
MARK_AS_ADVANCED(ILUPACK_INCLUDE_DIR)

