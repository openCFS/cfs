SET(CFXCUST_FOUND 0)

#-------------------------------------------------------------------------------
# Look for CFX Customizations header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("CFX Customizations"
  "${CFS_BINARY_DIR}/include/CFX_v100/cfx5ext.h"
  "${CFS_DEPS_ROOT}/cfx_custom/build_cfxcust.pl"
  "build_cfxcust.log")

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

SET(CFX_V100_INCLUDE_DIR "${CFS_BINARY_DIR}/include/CFX_v100")
SET(CFX_V110_INCLUDE_DIR "${CFS_BINARY_DIR}/include/CFX_v110")

MARK_AS_ADVANCED(CFX_V100_LIBIO)
MARK_AS_ADVANCED(CFX_V110_LIBIO)

SET(CFXCUST_FOUND 1)
