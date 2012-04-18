#-------------------------------------------------------------------------------
# SCPIP - - an efficient software tool for the solution of
# structural optimization problems.
#
# Project Homepage
# http://www.mathematik.uni-wuerzburg.de/~zillober
# http://www.mathematik.uni-wuerzburg.de/~zillober/pubs/manual30.pdf
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to scpip sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(scpip_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/scpip")
set(scpip_source  "${scpip_prefix}/src/scpip")
set(scpip_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${scpip_install}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
#  -DCMAKE_Fortran_FLAGS:FILEPATH=-v
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
)

#-------------------------------------------------------------------------------
# The scpip external project
#-------------------------------------------------------------------------------
ExternalProject_Add(scpip
  PREFIX "${scpip_prefix}"
  SOURCE_DIR "${scpip_source}"
  URL ${SCPIP_PATH}/${SCPIP_BZ2}
  URL_MD5 ${SCPIP_MD5}
  CMAKE_ARGS
    ${CMAKE_ARGS}
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/scpip/scpip-patch.cmake.in")
SET(PFN "${scpip_prefix}/scpip-patch.cmake")
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
ExternalProject_Add_Step(scpip custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${scpip_source}
)


SET(SCPIP_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of SCPIP libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(SCPIP_LIBRARY
  "${LD}/libDscpip.a;${LD}/libZscpip.a;${LD}/libblaslike.a;${LD}/libsparspak.a;${LD}/libamd.a"
  CACHE FILEPATH "SCPIP library.")

MARK_AS_ADVANCED(SCPIP_LIBRARY)
MARK_AS_ADVANCED(SCPIP_INCLUDE_DIR)

