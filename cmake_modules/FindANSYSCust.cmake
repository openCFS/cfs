SET(ANSCUST_FOUND 0)

#-------------------------------------------------------------------------------
# Look for ANSYS Customization header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("ANSYS Customizations"
  "${CFS_BINARY_DIR}/include/v110/ansys.h"
  "${CFS_DEPS_ROOT}/ansys_custom/build_anscust.pl"
  "build_anscust.log")

#-------------------------------------------------------------------------------
# Determine paths of MUPARSER libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")

IF(CFS_ARCH STREQUAL "I386")
  SET(ANSYS_V100_LIBBIN
    "${LD}/v100/linia32/libbin.so"
    CACHE FILEPATH "ANSYS binary interface library.")
  SET(ANSYS_V110_LIBBIN
    "${LD}/v110/linia32/libbin.so"
    CACHE FILEPATH "ANSYS binary interface library.")
ENDIF(CFS_ARCH STREQUAL "I386")

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
    "${LD}/v110/linem64t/libbin.so"
    CACHE FILEPATH "ANSYS binary interface library.")
ENDIF(CFS_SUBARCH STREQUAL "EM64T")

SET(ANSYS_V100_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v100")
SET(ANSYS_V110_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v110")

MARK_AS_ADVANCED(ANSYS_V100_LIBBIN)
MARK_AS_ADVANCED(ANSYS_V110_LIBBIN)

SET(ANSCUST_FOUND 1)
