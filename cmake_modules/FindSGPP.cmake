# - Find SGPP installation 
# This module determines where the SGPP libraries are.
#
#  SGPP_FOUND   = system has SGPP
#  SGPP_LIBRARY = path to the SGPP libraries
#
# Daniel Hübner (03/14)
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SET(SGPP_FOUND 0)

#-------------------------------------------------------------------------------
# Look for SGPP header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("SGpp"
  "${CFS_BINARY_DIR}/include/sgpp_base.hpp"
  "${CFS_DEPS_ROOT}/sgpp/build_sgpp.pl"
  "build_sgpp.log")

SET(SGPP_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of SGPP libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(SGPP_LIBRARY "${LD}/libsgppopt.a;${LD}/libumfpack.a;${LD}/libsuitesparseconfig.a;${LD}/libarmadillo.a;${LD}/libsgppbase.a;${LD}/libsgpppde.a;${LD}/libsgppsolver.a"
    CACHE FILEPATH "SGPP library.")

MARK_AS_ADVANCED(SGPP_INCLUDE_DIR)
MARK_AS_ADVANCED(SGPP_LIBRARY)

# Add a switch for the compiler.
# As the SGpp toolbox is not available for everyone, we have to exclude the corresponding source code from compiling.
IF(USE_SGPP)
  ADD_DEFINITIONS(-DUSE_SGPP)
ENDIF()

SET(SGPP_FOUND 1)
