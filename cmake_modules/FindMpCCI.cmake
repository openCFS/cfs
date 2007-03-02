# - Find glib2 installation 
# This module finds if glib2 is installed and determines where 
# the include files and libraries are. It also determines what the name of
# the library is. This code sets the following variables:
#  
#  MPCCI_FOUND     = system has glib2 
#  MPCCI_LIBRARIES = path to the glib2 libraries
#  MPCCI_CXX_FLAGS  = Compiler flags for glib2 
#  MPCCI_INCLUDE_DIR      = where to find "glib2.h"
#  MPCCI_DEFINITIONS      = extra defines
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

OPTION(USE_MPCCI_305 "Use the new (3.0.5) or the old (3.0.3) MpCCI" ON)

IF(WIN32_STYLE_FIND)

  ## ######################################################################
  ##
  ## Windows specific:
  ##
  ## candidates for root/base directory of Glib2
  ## should have subdirs include and lib containing source/glib2.h
  ## fix the root dir to avoid mixing of headers/libs from different
  ## versions/builds:
  
  SET (MPCCI_POSSIBLE_ROOT_PATHS
    $ENV{MPCCI_ROOT}
    "C:\\Dev\\glib2"
    )
  
  FIND_PATH(MPCCI_ROOT_DIR  source/glib2.h 
    ${MPCCI_POSSIBLE_ROOT_PATHS} )  
  
  
  ## find libs for combination of static/shared with release/debug
  ## be careful if you add something here, 
  SET (MPCCI_POSSIBLE_LIB_PATHS
    "${MPCCI_ROOT_DIR}/win/Release"
    "${MPCCI_ROOT_DIR}/win/Debug"
    ) 
  
  FIND_LIBRARY(MPCCI_STATIC_LIBRARY
    NAMES glib2
    PATHS 
    "${MPCCI_ROOT_DIR}/win/Release"
    ${MPCCI_POSSIBLE_LIB_PATHS}
    DOC "glib2 static release build library" ) 
  
  FIND_LIBRARY(MPCCI_STATIC_DEBUG_LIBRARY
    NAMES glib2
    PATHS 
    "${MPCCI_ROOT_DIR}/win/Debug"
    ${MPCCI_POSSIBLE_LIB_PATHS}
    DOC "glib2 static debug build library" )

  
  
  ##
  ## now we should have found all Glib2 libs available on the system.
  ## let the user decide which of the available onse to use.
  ## 
  
  SET(MPCCI_LIBRARIES
    debug ${MPCCI_STATIC_DEBUG_LIBRARY}   optimized ${MPCCI_STATIC_LIBRARY}
    )


  
  
  
  
  ## Find the include directories for Glib2
  ## add inc dir for general for "glib2.h"
  FIND_PATH(MPCCI_INCLUDE_DIR  glib2.h 
    "${MPCCI_ROOT_DIR}/source" )  
  
  MARK_AS_ADVANCED(
    MPCCI_ROOT_DIR
    MPCCI_INCLUDE_DIR
    MPCCI_STATIC_LIBRARY
    MPCCI_STATIC_DEBUG_LIBRARY
    )
  
  
ELSE(WIN32_STYLE_FIND)

  IF (UNIX_STYLE_FIND) 

    IF(USE_MPCCI_305)
      SET (MPICH_POSSIBLE_ROOT_PATHS
          "$ENV{CFSDEPS_ROOT}"
       	  "/opt/mpich/ch-p4"
          #          "/home/data/programs/MpCCI/MpCCI_3.0.5/AMD64/MPICH64/mpich-1.2.7p1"
          #	"/home/data/programs/mpich-1.2.7"
          #	"/opt/mpich/ch-p4mpd"
	)

      
      FIND_PATH(MPICH_ROOT_DIR include/mpi.h 
	${MPICH_POSSIBLE_ROOT_PATHS} )
      
      SET (MPICH_POSSIBLE_LIB_PATHS
	"${MPICH_ROOT_DIR}/lib64"
	"${MPICH_ROOT_DIR}/lib"
	) 
      
      FIND_LIBRARY(MPICH_LIBRARY
	NAMES mpich
	PATHS ${MPICH_POSSIBLE_LIB_PATHS}
	DOC "MPICH library" ) 
      
      SET(MPICH_LIBRARIES
	"${MPICH_LIBRARY}"
	) 
      
      #/opt/mpich/ch-p4/lib64/libmpich.a
#      MESSAGE("DBG found MPICH_LIBRARIES: ${MPICH_LIBRARIES}")
      
      STRING(REGEX REPLACE "/libmpich.*" #"[^/]*$"
        "" MPICH_LIBRARY_DIR
        ${MPICH_LIBRARIES})
      
#      MESSAGE("DBG found MPICH_LIBRARY_DIR: ${MPICH_LIBRARY_DIR}")
      
      
      SET(MPCCI_ROOT_DIR "/home/data/programs/MpCCI/MpCCI_3.0.5")

      #    MESSAGE("DBG found MPCCI_ROOTDIR: ${MPCCI_ROOT_DIR}")

      EXEC_PROGRAM("${MPCCI_ROOT_DIR}/bin/mpcci"
	ARGS arch -n
	OUTPUT_VARIABLE MPCCI_ARCH
	RETURN_VALUE RETVAL)

      IF(MPCCI_ARCH STREQUAL "linux_amd64")
	SET(MPCCI_LIBRARY -lmpccisdk-64)
      ENDIF(MPCCI_ARCH STREQUAL "linux_amd64")

      IF(MPCCI_ARCH STREQUAL "linux_em64t")
	SET(MPCCI_LIBRARY -lmpccisdk-64)
      ENDIF(MPCCI_ARCH STREQUAL "linux_em64t")

      IF(MPCCI_ARCH STREQUAL "linux_x86")
	SET(MPCCI_LIBRARY -lmpccisdk-32)
      ENDIF(MPCCI_ARCH STREQUAL "linux_x86")

      SET(MPCCI_LIBRARIES "-L${MPCCI_ROOT_DIR}/lib/${MPCCI_ARCH};${MPCCI_LIBRARY};-L${MPICH_LIBRARY_DIR};${MPICH_LIBRARY}")
      SET(MPCCI_INCLUDE_DIR "${MPCCI_ROOT_DIR}/include;${MPICH_ROOT_DIR}/include")
      SET(MPCCI_CXX_FLAGS "-fPIC -O3 -Wall")

      EXEC_PROGRAM("${MPCCI_ROOT_DIR}/bin/mpcci"
	ARGS info -release
	OUTPUT_VARIABLE MPCCI_RELEASE
	RETURN_VALUE RETVAL)

#      MESSAGE("DBG found MPCCI_RELEASE: ${MPCCI_RELEASE}")
    ELSE(USE_MPCCI_305)
      SET(MPCCI_ROOT_DIR "/home/data/programs/MpCCI/mpcci-3.0.3-sdk-linux-x86")
      SET(MPICH_ROOT_DIR "/home/data/programs/mpich-1.2.7")

      SET(MPCCI_LIBRARIES "-L${MPCCI_ROOT_DIR}/lib;-lcci;-L${MPICH_ROOT_DIR}/lib;-lmpich;-lg2c;-lstdc++")
      SET(MPCCI_INCLUDE_DIR "${MPCCI_ROOT_DIR}/include;${MPICH_ROOT_DIR}/include")
      SET(MPCCI_CXX_FLAGS "")

      EXEC_PROGRAM("${MPCCI_ROOT_DIR}/bin/ccirun"
	ARGS -version
	OUTPUT_VARIABLE MPCCI_RELEASE
	RETURN_VALUE RETVAL)

    ENDIF(USE_MPCCI_305)
    
    STRING(REGEX MATCH "3\\.0\\.[0-9]+"
      MPCCI_RELEASE
      ${MPCCI_RELEASE})

    STRING(REGEX REPLACE "\\." ""
      MPCCI_RELEASE
      ${MPCCI_RELEASE})

#    MESSAGE("DBG found MPCCI_RELEASE: ${MPCCI_RELEASE}")

#    MESSAGE("DBG found MPCCI_LIBRARIES: ${MPCCI_LIBRARIES}")
#    MESSAGE("DBG found MPCCI_INCLUDE_DIR: ${MPCCI_INCLUDE_DIR}")
#    MESSAGE("DBG found MPCCI_CXX_FLAGS: ${MPCCI_CXX_FLAGS}")
#    MESSAGE("DBG found MPCCI_LD_FLAGS: ${MPCCI_LD_FLAGS}")

  ELSE(UNIX_STYLE_FIND)
    MESSAGE(STATUS "FindMpCCI.cmake:  Platform unknown/unsupported by FindGlib2.cmake. It's neither WIN32 nor UNIX")
  ENDIF(UNIX_STYLE_FIND)
ENDIF(WIN32_STYLE_FIND)


IF(MPCCI_LIBRARIES)
  IF(MPCCI_INCLUDE_DIR OR MPCCI_CXX_FLAGS)
    ## found all we need.
    SET(MPCCI_FOUND 1)
    SET(MPCCI_INCLUDE_DIR "${MPCCI_INCLUDE_DIR}"
        CACHE INTERNAL "MPCCI_INCLUDE_DIR")
    SET(MPCCI_CXX_FLAGS "${MPCCI_CXX_FLAGS}"
        CACHE INTERNAL "MPCCI_CXX_FLAGS")
    SET(MPCCI_LIBRARIES "${MPCCI_LIBRARIES};-lnsl;-lm"
        CACHE INTERNAL "MPCCI_LIBRARIES")
    SET(MPCCI_RELEASE "${MPCCI_RELEASE}"
        CACHE INTERNAL "MPCCI_RELEASE")
    
  ENDIF(MPCCI_INCLUDE_DIR OR MPCCI_CXX_FLAGS)
ENDIF(MPCCI_LIBRARIES)

MESSAGE("DBG found MPCCI_FOUND: ${MPCCI_FOUND}")
