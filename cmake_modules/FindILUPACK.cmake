SET(ILUPACK_FOUND 0)

#-------------------------------------------------------------------------------
# Look for ILUPACK header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("AMD"
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libamd.a"
  "${CFS_DEPS_ROOT}/amd/build_amd.pl"
  "build_amd.log")

BUILD_EXTLIB("ILUPACK"
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libZilupack.a"
  "${CFS_DEPS_ROOT}/ilupack/build_ilupack.pl"
  "build_ilupack.log")

SET(ILUPACK_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of ILUPACK libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(ILUPACK_LIBRARY
  "${LD}/libDilupack.a;${LD}/libZilupack.a;${LD}/libblaslike.a;${LD}/libsparspak.a;${LD}/libamd.a"
  CACHE FILEPATH "ILUPACK library.")

MARK_AS_ADVANCED(ILUPACK_LIBRARY)
SET(ILUPACK_FOUND 1)
