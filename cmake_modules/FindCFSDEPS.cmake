# - Find CFSDEPS installation 
# This module finds if CFSDEPS is installed and determines where 
# the include files and libraries are. It also determines what the name of
# the library is. This code sets the following variables:
#  
#  CFSDEPS_FOUND     = system has CFSDEPS 
#  CFSDEPS_LIBRARIES = path to the CFSDEPS libraries
#  CFSDEPS_CXX_FLAGS  = Compiler flags for CFSDEPS 
#  CFSDEPS_INCLUDE_DIR      = where to find "CFSDEPS.h"
#  CFSDEPS_DEFINITIONS      = extra defines
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

  ## ######################################################################
  ##
  ## Windows specific:
  ##
  ## candidates for root/base directory of CFSDEPS
  ## should have subdirs include and lib containing source/CFSDEPS.h
  ## fix the root dir to avoid mixing of headers/libs from different
  ## versions/builds:
  
  SET (CFSDEPS_POSSIBLE_ROOT_PATHS
    $ENV{CFSDEPS_ROOT}
    "C:\\Dev\\CFSDEPS"
    )
  
  FIND_PATH(CFSDEPS_ROOT_DIR  source/CFSDEPS.h 
    ${CFSDEPS_POSSIBLE_ROOT_PATHS} )  
  
  
  ## find libs for combination of static/shared with release/debug
  ## be careful if you add something here, 
  SET (CFSDEPS_POSSIBLE_LIB_PATHS
    "${CFSDEPS_ROOT_DIR}/win/Release"
    "${CFSDEPS_ROOT_DIR}/win/Debug"
    ) 
  
  FIND_LIBRARY(CFSDEPS_STATIC_LIBRARY
    NAMES CFSDEPS
    PATHS 
    "${CFSDEPS_ROOT_DIR}/win/Release"
    ${CFSDEPS_POSSIBLE_LIB_PATHS}
    DOC "CFSDEPS static release build library" ) 
  
  FIND_LIBRARY(CFSDEPS_STATIC_DEBUG_LIBRARY
    NAMES CFSDEPS
    PATHS 
    "${CFSDEPS_ROOT_DIR}/win/Debug"
    ${CFSDEPS_POSSIBLE_LIB_PATHS}
    DOC "CFSDEPS static debug build library" )

  
  
  ##
  ## now we should have found all CFSDEPS libs available on the system.
  ## let the user decide which of the available onse to use.
  ## 
  
  SET(CFSDEPS_LIBRARIES
    debug ${CFSDEPS_STATIC_DEBUG_LIBRARY}   optimized ${CFSDEPS_STATIC_LIBRARY}
    )


  
  
  
  
  ## Find the include directories for CFSDEPS
  ## add inc dir for general for "CFSDEPS.h"
  FIND_PATH(CFSDEPS_INCLUDE_DIR  CFSDEPS.h 
    "${CFSDEPS_ROOT_DIR}/source" )  
  
  MARK_AS_ADVANCED(
    CFSDEPS_ROOT_DIR
    CFSDEPS_INCLUDE_DIR
    CFSDEPS_STATIC_LIBRARY
    CFSDEPS_STATIC_DEBUG_LIBRARY
    )
  
  
ELSE(WIN32_STYLE_FIND)

  IF (UNIX_STYLE_FIND) 

    SET (CFSDEPS_POSSIBLE_ROOT_PATHS
      "$ENV{CFSDEPS_ROOT}"
      "/home/data/libraries/CFSDEPS/${CFS_ARCH_STR}"
      "/opt/lse"
      )

#    MESSAGE("DBG CFSDEPS_POSSIBLE_ROOT_PATHS: ${CFSDEPS_POSSIBLE_ROOT_PATHS}")
    
    FIND_PATH(CFSDEPS_ROOT_DIR  include/hdf5.h 
      ${CFSDEPS_POSSIBLE_ROOT_PATHS} )  
    MARK_AS_ADVANCED(CFSDEPS_ROOT_DIR)
    
    SET (CFSDEPS_POSSIBLE_LIB_PATHS
      "${CFSDEPS_ROOT_DIR}/lib64"
      "${CFSDEPS_ROOT_DIR}/lib"
      ) 
    
    FIND_LIBRARY(TEST_HDF5_LIBRARY
      NAMES hdf5
      PATHS ${CFSDEPS_POSSIBLE_LIB_PATHS}
      DOC "HDF5 library" ) 
    MARK_AS_ADVANCED(TEST_HDF5_LIBRARY)

    STRING(REGEX REPLACE "/libhdf5.*" #"[^/]*$"
      "" CFSDEPS_LIBRARY_DIR
      ${TEST_HDF5_LIBRARY})

    FIND_PATH(CFSDEPS_INCLUDE_DIR hdf5.h 
      "${CFSDEPS_ROOT_DIR}/include" ) 
    MARK_AS_ADVANCED(CFSDEPS_INCLUDE_DIR)


  ELSE(UNIX_STYLE_FIND)
    MESSAGE(FATAL_ERROR "FindCFSDEPS.cmake:  Platform unknown/unsupported by FindCFSDEPS.cmake. It's neither WIN32 nor UNIX")
  ENDIF(UNIX_STYLE_FIND)
ENDIF(WIN32_STYLE_FIND)


IF(CFSDEPS_INCLUDE_DIR OR CFSDEPS_CXX_FLAGS)
  ## found all we need.
  SET(CFSDEPS_FOUND 1)
  SET(CFSDEPS_INCLUDE_DIR "${CFSDEPS_INCLUDE_DIR}"
    CACHE PATH INTERNAL "CFSDEPS_INCLUDE_DIR")
  
  SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS}"
    CACHE INTERNAL "CFSDEPS_CXX_FLAGS")
  
  SET(CFSDEPS_LIBRARIES "${CFSDEPS_LIBRARY_DIR}"
    CACHE INTERNAL "CFSDEPS_LIBRARY_DIR")
  
#    MESSAGE("CFSDEPS_INCLUDE_DIR: ${CFSDEPS_INCLUDE_DIR}")
#    MESSAGE("CFSDEPS_LIBRARY_DIR: ${CFSDEPS_LIBRARY_DIR}")
#    MESSAGE("DBG found CFSDEPS_FOUND: ${CFSDEPS_FOUND}")

ENDIF(CFSDEPS_INCLUDE_DIR OR CFSDEPS_CXX_FLAGS)
