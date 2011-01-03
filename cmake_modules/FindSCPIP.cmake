# - Find SCPIP installation 
# This module finds if SCPIP is installed and determines where 
# the libraries are.
#
#  SCPIP_FOUND     = system has SCPIP 
#  SCPIP_LIBRARY = path to the SCPIP libraries
#
# Fabian Wein (11/07)
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that one needs an academic use license from 
# PD. Christian Zillober (Uni-Würzburg) before using SCPIP!!

SET(SCPIP_FOUND 0)

IF(EXISTS "${CFS_DEPS_ROOT}/scpip/build_scpip.pl")
  BUILD_EXTLIB("SCPIP"
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libscpip.a"
    "${CFS_DEPS_ROOT}/scpip/build_scpip.pl"
    "build_scpip.log")

  #-----------------------------------------------------------------------------
  # Determine paths of SCPIP libraries.
  #-----------------------------------------------------------------------------
  SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
  SET(SCPIP_LIBRARY
    "${LD}/libscpip.a"
    CACHE FILEPATH "SCPIP library.")

  MARK_AS_ADVANCED(SCPIP_LIBRARY)

  SET(SCPIP_FOUND 1)
ELSE(EXISTS "${CFS_DEPS_ROOT}/scpip/build_scpip.pl")
  SET(USE_SCPIP OFF CACHE BOOL "Enable SCPIP" FORCE)
  MESSAGE(STATUS "SCPIP not found in CFSDEPS. Forcing it off!")
ENDIF(EXISTS "${CFS_DEPS_ROOT}/scpip/build_scpip.pl")
