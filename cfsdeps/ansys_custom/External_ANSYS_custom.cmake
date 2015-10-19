#-------------------------------------------------------------------------------
# ANSYS I/O libraries for reading and writing .rst files.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to ansys_custom sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(ansys_custom_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ansys_custom")
set(ansys_custom_source  "${ansys_custom_prefix}/src/ansys_custom")
set(ansys_custom_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set names of build file and template file.
#-------------------------------------------------------------------------------
SET(BUILD_TEMPL
  "${CFS_SOURCE_DIR}/cfsdeps/ansys_custom/ansys_custom-build.cmake.in")
SET(BUILD "${ansys_custom_prefix}/ansys_custom-build.cmake")
CONFIGURE_FILE("${BUILD_TEMPL}" "${BUILD}" @ONLY) 

#-------------------------------------------------------------------------------
# Set names of install file and template file.
#-------------------------------------------------------------------------------
SET(INST_TEMPL
  "${CFS_SOURCE_DIR}/cfsdeps/ansys_custom/ansys_custom-install.cmake.in")
SET(INST "${ansys_custom_prefix}/ansys_custom-install.cmake")
CONFIGURE_FILE("${INST_TEMPL}" "${INST}" @ONLY) 


#-------------------------------------------------------------------------------
# Configure target for downloading ANSYS customization archives using CMake
# ExternalData mechanism.
#-------------------------------------------------------------------------------
# Add standard remote object stores to user's configuration.
SET(ExternalData_URL_TEMPLATES
  "${CFS_DS_WEBDAV}/cfsdeps/sources/ansys_custom/%(algo)/%(hash)"
)

# Set standard local object stores.
SET(ExternalData_OBJECT_STORES
  "${CFS_DEPS_CACHE_DIR}/sources/ansys_custom"
)

SET(ARCHIVES
  "v100_include"
  "v110_include"
  "v120_include"
  "v121_include"
  "v130_include"
  "v140_include"

  "v100_linop64"
  "v110_linop64"

  "v100_linem64t"
  "v110_linem64t"

  "v120_linx64"
  "v121_linem64t"
  "v130_linx64"
  "v140_linx64"

  "v110_linia64"

  "v100_linia32"
  "v110_linia32"
  "v121_linia32"
  
  "v110_winx64"
  "v120_winx64"
  "v121_winx64"
  "v130_winx64"
  "v140_winx64"
  
  "v100_win32"
  "v110_win32"
  "v120_win32"
  "v121_win32"
  "v130_win32"
  "v140_win32"
)

foreach(archive ${ARCHIVES})
  set(ANSYS_CUSTOM_EXTERNAL_DATA
    "${ANSYS_CUSTOM_EXTERNAL_DATA}
DATA{cfsdeps/ansys_custom/ansys_custom_${archive}.tar.bz2}")

  # Give a hint about downloading the source archive to the developer.
  FILE(READ "cfsdeps/ansys_custom/ansys_custom_${archive}.tar.bz2.md5" ANSCUST_HASH)
  STRING(STRIP ${ANSCUST_HASH} ANSCUST_HASH)

  IF(NOT EXISTS "${ExternalData_OBJECT_STORES}/MD5/${ANSCUST_HASH}")
    SET(MSG "Please download the file ")
    SET(MSG "${MSG}'${CFS_DS_WEBDAV}/cfsdeps/sources/ansys_custom/MD5/${ANSCUST_HASH}'")
    SET(MSG "${MSG} to '${ExternalData_OBJECT_STORES}/MD5/${ANSCUST_HASH}'.")
    
    colormsg(HIYELLOW "${MSG}")
  ENDIF()

endforeach(archive)

# Expand all arguments as a single string to preserve escaped semicolons.
ExternalData_expand_arguments(ansys_custom_external_data
  targetArgs "${ANSYS_CUSTOM_EXTERNAL_DATA}")

# Add a build target to populate the real data.
ExternalData_Add_Target(ansys_custom_external_data)

# There is no ANSYS_CUSTOM_VER yet! 
PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "ansys_custom" "no_ver") 
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${ansys_custom_prefix}")

SET(ZIPFROMCACHE "${ansys_custom_prefix}/ansys_custom-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${ansys_custom_prefix}/ansys_custom-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The ansys_custom external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(ansys_custom
    PREFIX "${ansys_custom_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(ansys_custom
    DEPENDS ansys_custom_external_data
    PREFIX "${ansys_custom_prefix}"
    SOURCE_DIR "${ansys_custom_source}"
    BUILD_IN_SOURCE 1
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND  ${CMAKE_COMMAND} -P ${BUILD}
    INSTALL_COMMAND  ${CMAKE_COMMAND} -P ${INST}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(ansys_custom cfsdeps_zipToCache
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
SET(CFSDEPS ${CFSDEPS} ansys_custom)

#-------------------------------------------------------------------------------
# Determine paths of ANSYS Customization libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")

IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")

IF(CFS_ARCH STREQUAL "I386")
  SET(ANSYS_V100_LIBBIN
    "${LD}/v100/linia32/libbin.so"
    CACHE FILEPATH "ANSYS binary interface library.")
  SET(ANSYS_V110_LIBBIN
    "${LD}/v110/linia32/libbin.so"
    CACHE FILEPATH "ANSYS binary interface library.")
ENDIF(CFS_ARCH STREQUAL "I386")

IF(CFS_ARCH STREQUAL "X86_64")
  IF(CFS_SUBARCH STREQUAL "OPTERON")
    SET(ANSYS_V100_LIBBIN
      "${LD}/v100/linop64/libbin.so"
      CACHE FILEPATH "ANSYS binary interface library.")
    SET(ANSYS_V110_LIBBIN
      "${LD}/v110/linop64/libbin.so"
      CACHE FILEPATH "ANSYS binary interface library.")
  ENDIF(CFS_SUBARCH STREQUAL "OPTERON")

  IF(CFS_SUBARCH STREQUAL "EM64T")
    SET(ANSYS_V100_LIBBIN
      "${LD}/v100/linem64t/libbin.so"
      CACHE FILEPATH "ANSYS binary interface library.")
    SET(ANSYS_V110_LIBBIN
      "${LD}/v110/linem64t/libbin.so;-lguide;-lpthread"
      CACHE FILEPATH "ANSYS binary interface library.")
  ENDIF(CFS_SUBARCH STREQUAL "EM64T")
  
  SET(ANSYS_V120_LIBBIN
    "${LD}/v120/libbin.so;-lpthread"
    CACHE FILEPATH "ANSYS binary interface library.")

  SET(ANSYS_V121_LIBBIN
    "${LD}/v121/linem64t/libbin.so;-lguide;-lpthread"
    CACHE FILEPATH "ANSYS binary interface library.")

  SET(ANSYS_V130_LIBBIN
    "${LD}/v130/libbin.so;-lpthread"
    CACHE FILEPATH "ANSYS binary interface library.")
  
  SET(ANSYS_V140_LIBBIN
    "${LD}/v140/libbin.so;-lpthread"
    CACHE FILEPATH "ANSYS binary interface library.")

ENDIF(CFS_ARCH STREQUAL "X86_64")
ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")

IF(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  IF(NOT MINGW)
    IF(CFS_ARCH STREQUAL "X86_64")
      SET(ANSYS_V110_LIBBIN
        "${LD}/v110/winx64/binlib.lib"
        CACHE FILEPATH "ANSYS binary interface library.")

      SET(ANSYS_V120_LIBBIN
        "${LD}/v120/winx64/binlib.lib"
        CACHE FILEPATH "ANSYS binary interface library.")
    
      SET(ANSYS_V121_LIBBIN
        "${LD}/v121/winx64/libbin.lib"
        CACHE FILEPATH "ANSYS binary interface library.")
    
      SET(ANSYS_V130_LIBBIN
        "${LD}/v130/winx64/libbin.lib"
        CACHE FILEPATH "ANSYS binary interface library.")
    
      SET(ANSYS_V140_LIBBIN
        "${LD}/v140/winx64/libbin.lib"
        CACHE FILEPATH "ANSYS binary interface library.")
    ENDIF(CFS_ARCH STREQUAL "X86_64")
  ELSE(NOT MINGW)
    IF(CFS_TARGET_ARCH STREQUAL "X86_64")
      SET(ANSYS_V110_LIBBIN
        "${LD}/v110/winx64/libmingw_binlib.dll.a"
        CACHE FILEPATH "ANSYS binary interface library.")

      SET(ANSYS_V120_LIBBIN
        "${LD}/v120/winx64/libmingw_binlib.dll.a"
        CACHE FILEPATH "ANSYS binary interface library.")
    
      SET(ANSYS_V121_LIBBIN
        "${LD}/v121/winx64/libmingw_binlib.dll.a"
        CACHE FILEPATH "ANSYS binary interface library.")
    
      SET(ANSYS_V130_LIBBIN
        "${LD}/v130/winx64/libmingw_binlib.dll.a"
        CACHE FILEPATH "ANSYS binary interface library.")
    
      SET(ANSYS_V140_LIBBIN
        "${LD}/v140/winx64/libmingw_binlib.dll.a"
        CACHE FILEPATH "ANSYS binary interface library.")
    ENDIF(CFS_TARGET_ARCH STREQUAL "X86_64")
  ENDIF(NOT MINGW)
ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Windows")

SET(ANSYS_V100_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v100")
SET(ANSYS_V110_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v110")
SET(ANSYS_V120_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v120")
SET(ANSYS_V121_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v121")
SET(ANSYS_V130_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v130")
SET(ANSYS_V140_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v140")

MARK_AS_ADVANCED(ANSYS_V100_LIBBIN)
MARK_AS_ADVANCED(ANSYS_V110_LIBBIN)
MARK_AS_ADVANCED(ANSYS_V120_LIBBIN)
MARK_AS_ADVANCED(ANSYS_V121_LIBBIN)
MARK_AS_ADVANCED(ANSYS_V130_LIBBIN)
MARK_AS_ADVANCED(ANSYS_V140_LIBBIN)

