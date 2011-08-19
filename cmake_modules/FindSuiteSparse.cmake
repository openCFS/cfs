SET(CHOLMOD_FOUND 0)

#-------------------------------------------------------------------------------
# Look for CholMod header.
#-------------------------------------------------------------------------------
CONFIGURE_FILE("${CFS_DEPS_ROOT}/suitesparse/build_suitesparse_vars.pl.in"
  "${CFS_TEMP_DIR}/build_suitesparse_vars.pl"
  @ONLY )

BUILD_EXTLIB("AMD/CHOLMOD from SuiteSparse"
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libcholmod.a"
  "${CFS_DEPS_ROOT}/suitesparse/build_suitesparse.pl"
  "build_suitesparse.log")

SET(CHOLMOD_INCLUDE_DIR "${CFSDEPS_INCLUDE_DIR}")

#-------------------------------------------------------------------------------
# Determine paths of CholMod libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(CHOLMOD_LIBRARY
  "${LD}/libcholmod.a;${LD}/libamd.a;${LD}/libcolamd.a;${LD}/libcamd.a;${LD}/libccolamd.a"
  CACHE FILEPATH "CholMod library.")

MARK_AS_ADVANCED(CHOLMOD_LIBRARY)
SET(CHOLMOD_FOUND 1)
