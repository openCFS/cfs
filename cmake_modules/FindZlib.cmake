#-------------------------------------------------------------------------------
# Find Zlib.
#-------------------------------------------------------------------------------
# FIND_PACKAGE(ZLIB)

MARK_AS_ADVANCED(ZLIB_LIBRARY)
MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# Look for zlib header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("zlib"
  "${CFS_BINARY_DIR}/include/zlib.h"
  "${CFS_DEPS_ROOT}/zlib/build_zlib.pl"
  "build_zlib.log")

SET(ZLIB_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libz.a CACHE FILEPATH "zlib library")
SET(ZLIB_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "zlib include directory")

MARK_AS_ADVANCED(ZLIB_LIBRARY)
MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR)

