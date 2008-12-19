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


SET(UNIX_STYLE_FIND 1)


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
