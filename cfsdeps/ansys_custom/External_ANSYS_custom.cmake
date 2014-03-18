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
# The ansys_custom external project
#-------------------------------------------------------------------------------
ExternalProject_Add(ansys_custom
  PREFIX "${ansys_custom_prefix}"
  SOURCE_DIR "${ansys_custom_source}"
  BUILD_IN_SOURCE 1
  DOWNLOAD_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND  ${CMAKE_COMMAND} -P ${BUILD}
  INSTALL_COMMAND  ${CMAKE_COMMAND} -P ${INST}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  ansys_custom
)

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

