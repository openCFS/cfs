SET(GIDPOST_FOUND 0)

#-------------------------------------------------------------------------------
# Determine path of GiDpost library.
#-------------------------------------------------------------------------------
FIND_LIBRARY(GIDPOST_LIBRARY
  NAMES gidpost
  PATHS ${CFSDEPS_LIBRARY_DIR}
  )

#-------------------------------------------------------------------------------
# Mark path of GIDPOST library as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(GIDPOST_LIBRARY)

FIND_PATH(GIDPOST_INCLUDE_DIR gidpost.h 
  ${CFSDEPS_INCLUDE_DIR} 
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH)

MARK_AS_ADVANCED(GIDPOST_INCLUDE_DIR)


IF(GIDPOST_LIBRARY AND GIDPOST_INCLUDE_DIR)
  SET(GIDPOST_FOUND 1)
ENDIF(GIDPOST_LIBRARY AND GIDPOST_INCLUDE_DIR)

IF(GIDPOST_FOUND)
  #-----------------------------------------------------------------------------
  # Determine version of GiDpost by reading its version header
  #-----------------------------------------------------------------------------
  FILE(READ "${CFSDEPS_INCLUDE_DIR}/gidpost.h" GIDPOST_HEADER)

  STRING(REGEX MATCH " gidpost [0-9]+\\.[0-9]+"
    CFS_GIDPOST_VERSION
    "${GIDPOST_HEADER}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+"
        CFS_GIDPOST_VERSION "${CFS_GIDPOST_VERSION}")

#  MESSAGE("CFS_GIDPOST_VERSION ${CFS_GIDPOST_VERSION}")

ENDIF(GIDPOST_FOUND)

MARK_AS_ADVANCED(GIDPOST_LIBRARY)

IF(NOT GIDPOST_FOUND)
  MESSAGE("Warning: GIDPOST could not be found! Please specify proper paths.")
ENDIF(NOT GIDPOST_FOUND)
