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
# Configure target for downloading SCPIP sources using CMake ExternalData
# mechanism.
#-------------------------------------------------------------------------------
# Add standard remote object stores to user's configuration.
list(APPEND ExternalData_URL_TEMPLATES
  "${CFS_DS_WEBDAV}/cfsdeps/sources/scpip/%(algo)/%(hash)"
  )

# Set standard local object stores.
SET(ExternalData_OBJECT_STORES
  "${CFS_DEPS_CACHE_DIR}/sources/scpip"
)

# Give a hint about downloading the source archive to the developer.
FILE(READ "cfsdeps/scpip/scpip.tar.bz2.md5" SCPIP_HASH)
STRING(STRIP ${SCPIP_HASH} SCPIP_HASH)

IF(NOT EXISTS "${ExternalData_OBJECT_STORES}/MD5/${SCPIP_HASH}")
  SET(MSG "Please download the file ")
  SET(MSG "${MSG}'${CFS_DS_WEBDAV}/cfsdeps/sources/scpip/MD5/${SCPIP_HASH}'")
  SET(MSG "${MSG} to '${ExternalData_OBJECT_STORES}/MD5/${SCPIP_HASH}'.")

  colormsg(HIYELLOW "${MSG}")
ENDIF()

set(SCPIP_EXTERNAL_DATA "DATA{cfsdeps/scpip/scpip.tar.bz2}")

# Expand all arguments as a single string to preserve escaped semicolons.
ExternalData_expand_arguments(scpip_external_data
  targetArgs "${SCPIP_EXTERNAL_DATA}")

# Add a build target to populate the real data.
ExternalData_Add_Target(scpip_external_data)

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${scpip_install}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:FILEPATH=-w
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
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/scpip/scpip-patch.cmake.in")
SET(PFN "${scpip_prefix}/scpip-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# The scpip external project
#-------------------------------------------------------------------------------
ExternalProject_Add(scpip
  DEPENDS scpip_external_data
  PREFIX "${scpip_prefix}"
  SOURCE_DIR "${scpip_source}"
  URL ${SCPIP_PATH}/${SCPIP_BZ2}
  URL_MD5 ${SCPIP_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CMAKE_ARGS
    ${CMAKE_ARGS}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  scpip
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

