# - Find KNITRO installation
# 
# AUTHOR
# Fabian Wein

SET(KNITRO_FOUND 0)

#-------------------------------------------------------------------------------
# Look for KNITRO header.
#-------------------------------------------------------------------------------
#BUILD_EXTLIB("KNITRO"
#  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libkniro700.a"
#  "${CFS_DEPS_ROOT}/ipopt/build_ipopt.pl"
#  "build_ipopt.log")

SET(KNITRO_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of KNITRO libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(KNITRO_LIBRARY
  "${LD}/libknitro700.a"
  CACHE FILEPATH "KNITRO library.")

MARK_AS_ADVANCED(KNITRO_LIBRARY)
SET(KNITRO_FOUND 1)

