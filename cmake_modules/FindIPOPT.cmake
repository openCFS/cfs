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
# Simon Triebenbacher simon@ibtriebenbacher.de (07/2006)


IF(WIN32)
  SET(WIN32_STYLE_FIND 1)
ENDIF(WIN32)
IF(MINGW)
  SET(WIN32_STYLE_FIND 0)
  SET(UNIX_STYLE_FIND 1)
ENDIF(MINGW)
IF(UNIX)
  SET(UNIX_STYLE_FIND 1)
ENDIF(UNIX)


IF(WIN32_STYLE_FIND)
# We do not use WIndows at all!
ELSE(WIN32_STYLE_FIND)

  IF (UNIX_STYLE_FIND) 

  SET (IPOPT_POSSIBLE_ROOT_PATHS
    $ENV{IPOPT_ROOT}
    "/space/fwein/packages/ipopt_3.3_stable/build"
    )
  
  FIND_PATH(IPOPT_ROOT_DIR include/config_ipopt.h
    ${IPOPT_POSSIBLE_ROOT_PATHS} )  

    
    SET (IPOPT_POSSIBLE_LIB_PATHS
      "${IPOPT_ROOT_DIR}/lib64"
      "${IPOPT_ROOT_DIR}/lib"
      ) 
    
    FIND_LIBRARY(IPOPT_LIBRARY
      NAMES ipopt
      PATHS ${IPOPT_POSSIBLE_LIB_PATHS}
      DOC "IPOPT library" ) 

    STRING(REGEX REPLACE "/libipopt.*" #"[^/]*$"
      "" IPOPT_LINK_DIR
      ${IPOPT_LIBRARY})

  ELSE(UNIX_STYLE_FIND)
    MESSAGE(FATAL_ERROR "FindIPOPT.cmake:  Platform unknown/unsupported by FindIPOPT.cmake. It's neither WIN32 nor UNIX")
  ENDIF(UNIX_STYLE_FIND)
ENDIF(WIN32_STYLE_FIND)


IF(IPOPT_ROOT_DIR)
  IF(IPOPT_LIBRARY)
  ## found all we need.
    SET(IPOPT_FOUND 1)

     SET(IPOPT_INCLUDE_DIR "${IPOPT_ROOT_DIR}/include"
       CACHE PATH "IPOPT_INCLUDE_DIR")

     SET(IPOPT_LIBRARY "${IPOPT_LIBRARY}"
       CACHE PATH "IPOPT_LIBRARY")

     SET(IPOPT_LINK_DIR "${IPOPT_LINK_DIR}"
       CACHE PATH "IPOPT_LINK_DIR")

     MARK_AS_ADVANCED(
       IPOPT_INCLUDE_DIR
       IPOPT_LIBRARY
       IPOPT_ROOT_DIR
     )
  ENDIF(IPOPT_LIBRARY)
ENDIF(IPOPT_ROOT_DIR)
