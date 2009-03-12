# - Find IPOPT installation 
# This module finds if IPOPT is installed and determines where 
# the include files and libraries are. It also determines what the name of
# the library is. This code sets the following variables:
#  
#  IPOPT_FOUND     = system has IPOPT 
#  IPOPT_LIBRARY = path to the IPOPT libraries
#  IPOPT_LINK_DIR = path to the IPOPT libraries
#  IPOPT_INCLUDE_DIR      = where to find "IPOPT.h"
#  IPOPT_ROOT_DIR      = where to find IPOPT
# 
# DEPRECATED
#
# OPTIONS 
# 
# USAGE 
# 
# NOTES
#
# AUTHOR
# Simon Triebenbacher simon.triebenbacher@uni-klu.ac.at (02/2009)

SET(IPOPT_FOUND 0)

#-------------------------------------------------------------------------------
# Look for IPOPT header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("IPOPT"
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libipopt.a"
  "${CFS_DEPS_ROOT}/ipopt/build_ipopt.pl"
  "build_ipopt.log")

SET(IPOPT_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of IPOPT libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(IPOPT_LIBRARY
  "${LD}/libipopt.a"
  CACHE FILEPATH "IPOPT library.")

MARK_AS_ADVANCED(IPOPT_LIBRARY)
SET(IPOPT_FOUND 1)

