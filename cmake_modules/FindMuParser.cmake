SET(MUPARSER_FOUND 0)

#-------------------------------------------------------------------------------
# Look for MUPARSER header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("MuParser"
  "${CFS_BINARY_DIR}/include/muParser.h"
  "${CFS_DEPS_ROOT}/muparser/build_muparser.pl"
  "build_muparser.log")

SET(MUPARSER_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of MUPARSER libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(MUPARSER_LIBRARY
  "${LD}/libmuparser.a"
  CACHE FILEPATH "MUPARSER library.")

MARK_AS_ADVANCED(MUPARSER_LIBRARY)
SET(MUPARSER_FOUND 1)
