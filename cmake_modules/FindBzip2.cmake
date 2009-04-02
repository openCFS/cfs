#-------------------------------------------------------------------------------
# Find bzip2.
#-------------------------------------------------------------------------------

MARK_AS_ADVANCED(BZIP2_LIBRARY)
MARK_AS_ADVANCED(BZIP2_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# Look for bzip2 header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("bzip2"
  "${CFS_BINARY_DIR}/include/bzlib.h"
  "${CFS_DEPS_ROOT}/bzip2/build_bzip2.pl"
  "build_bzip2.log")

SET(BZIP2_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libbz2.a CACHE FILEPATH "zlib library")
SET(BZIP2_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "bzip2 include directory")

MARK_AS_ADVANCED(BZIP2_LIBRARY)
MARK_AS_ADVANCED(BZIP2_INCLUDE_DIR)

