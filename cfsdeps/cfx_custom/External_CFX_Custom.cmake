#-------------------------------------------------------------------------------
# ANSYS CFX I/O libraries for reading .def, .res and .trn files.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to cfxio sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(cfxio_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cfxio")
set(cfxio_install  "${CMAKE_CURRENT_BINARY_DIR}")
set(cfxio_source  "${cfxio_prefix}/src/cfxio")

SET(BUILDCMD_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cfx_custom/cfxio-build.cmake.in")
SET(BUILDCMD "${cfxio_prefix}/cfxio-build.cmake")
CONFIGURE_FILE("${BUILDCMD_TEMPL}" "${BUILDCMD}" @ONLY) 

SET(INSTCMD_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cfx_custom/cfxio-install.cmake.in")
SET(INSTCMD "${cfxio_prefix}/cfxio-install.cmake")
CONFIGURE_FILE("${INSTCMD_TEMPL}" "${INSTCMD}" @ONLY) 

#-------------------------------------------------------------------------------
# Configure target for downloading CFX I/O archives using CMake
# ExternalData mechanism.
#-------------------------------------------------------------------------------
# Add standard remote object stores to user's configuration.
list(APPEND ExternalData_URL_TEMPLATES
  "${CFS_DS_WEBDAV}/cfsdeps/sources/cfx_custom/%(algo)/%(hash)"
  )

# Set standard local object stores.
SET(ExternalData_OBJECT_STORES
  "${CFS_DEPS_CACHE_DIR}/sources/cfx_custom"
)

SET(ARCHIVES
  "v100_include"
  "v110_include"
  "v120_include"
  "v121_include"
  "v130_include"
  "v140_include"

  "v100_linux-amd64"
  "v110_linux-amd64"
  "v120_linux-amd64"
  "v121_linux-amd64"
  "v130_linux-amd64"
  "v140_linux-amd64"
  
  "v110_linux-ia64"
  
  "v100_linux"
  "v110_linux"
  )

foreach(archive ${ARCHIVES})
  set(CFX_CUSTOM_EXTERNAL_DATA
    "${CFX_CUSTOM_EXTERNAL_DATA}
DATA{cfsdeps/cfx_custom/cfx_custom_${archive}.tar.bz2}")

  # Give a hint about downloading the source archive to the developer.
  FILE(READ "cfsdeps/cfx_custom/cfx_custom_${archive}.tar.bz2.md5" CFXCUST_HASH)
  STRING(STRIP ${CFXCUST_HASH} CFXCUST_HASH)

  IF(NOT EXISTS "${ExternalData_OBJECT_STORES}/MD5/${CFXCUST_HASH}")
    SET(MSG "Please download the file ")
    SET(MSG "${MSG}'${CFS_DS_WEBDAV}/cfsdeps/sources/cfx_custom/MD5/${CFXCUST_HASH}'")
    SET(MSG "${MSG} to '${ExternalData_OBJECT_STORES}/MD5/${CFXCUST_HASH}'.")
    
    colormsg(HIYELLOW "${MSG}")
  ENDIF()

endforeach(archive)

# Expand all arguments as a single string to preserve escaped semicolons.
ExternalData_expand_arguments(cfx_custom_external_data
  targetArgs "${CFX_CUSTOM_EXTERNAL_DATA}")

# Add a build target to populate the real data.
ExternalData_Add_Target(cfx_custom_external_data)

PRECOMPILED_ZIP_CXX(PRECOMPILED_PCKG_FILE "cfxio" "${CFXIO_VER}")

# This should be either PREFIX_DIR/src (install manifest is used for zipping)
# or PREFIX_DIR/install (install directory will be zipped)
SET(TMP_DIR "${cfxio_prefix}/src")

SET(ZIPFROMCACHE "${cfxio_prefix}/cfxio-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${cfxio_prefix}/cfxio-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The cfxio external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(cfxio
    PREFIX "${cfxio_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  ExternalProject_Add(cfxio
    DEPENDS cfx_custom_external_data
    PREFIX "${cfxio_prefix}"
    SOURCE_DIR "${cfxio_source}"
    BUILD_IN_SOURCE 1
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND  ${CMAKE_COMMAND} -P ${BUILDCMD}
    INSTALL_COMMAND  ${CMAKE_COMMAND} -P ${INSTCMD}
  )

  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(cfxio cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  cfxio
)

#-------------------------------------------------------------------------------
# Determine paths of CFX libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")

IF(CFS_ARCH STREQUAL "IA64" AND
   CFS_CXX_COMPILER_NAME STREQUAL "GCC")
  SET(CFX_IRC_LIB "${LD}/CFX_v110/libirc.a")
ENDIF(CFS_ARCH STREQUAL "IA64" AND
      CFS_CXX_COMPILER_NAME STREQUAL "GCC")

SET(CFX_V100_LIBIO
  "${LD}/CFX_v100/libio.a"
  CACHE FILEPATH "CFX binary interface library.")
SET(CFX_V110_LIBIO
  "${LD}/CFX_v110/libio.a;${CFX_IRC_LIB}"
  CACHE FILEPATH "CFX binary interface library.")
SET(CFX_V120_LIBIO
  "${LD}/CFX_v120/libio.a"
  CACHE FILEPATH "CFX binary interface library.")
SET(CFX_V121_LIBIO
  "${LD}/CFX_v121/libio.a"
  CACHE FILEPATH "CFX binary interface library.")
SET(CFX_V130_LIBIO
  "${LD}/CFX_v130/libio.a"
  CACHE FILEPATH "CFX binary interface library.")
SET(CFX_V140_LIBIO
  "${LD}/CFX_v140/libio.a"
  CACHE FILEPATH "CFX binary interface library.")


SET(CFX_V100_INCLUDE_DIR "${CFS_BINARY_DIR}/include/CFX_v100")
SET(CFX_V110_INCLUDE_DIR "${CFS_BINARY_DIR}/include/CFX_v110")
SET(CFX_V120_INCLUDE_DIR "${CFS_BINARY_DIR}/include/CFX_v120")
SET(CFX_V121_INCLUDE_DIR "${CFS_BINARY_DIR}/include/CFX_v121")
SET(CFX_V130_INCLUDE_DIR "${CFS_BINARY_DIR}/include/CFX_v130")
SET(CFX_V140_INCLUDE_DIR "${CFS_BINARY_DIR}/include/CFX_v140")

IF(CFS_ARCH STREQUAL "X86_64")
  SET(CFX_INCLUDE_DIR "${CFX_V140_INCLUDE_DIR}")
  SET(CFX_LIBIO "${CFX_V140_LIBIO}")
ELSE(CFS_ARCH STREQUAL "X86_64")
  SET(CFX_INCLUDE_DIR "${CFX_V110_INCLUDE_DIR}")
  SET(CFX_LIBIO "${CFX_V110_LIBIO}")
ENDIF(CFS_ARCH STREQUAL "X86_64")

MARK_AS_ADVANCED(CFX_V100_LIBIO)
MARK_AS_ADVANCED(CFX_V110_LIBIO)
MARK_AS_ADVANCED(CFX_V120_LIBIO)
MARK_AS_ADVANCED(CFX_V121_LIBIO)
MARK_AS_ADVANCED(CFX_V130_LIBIO)
MARK_AS_ADVANCED(CFX_V140_LIBIO)
