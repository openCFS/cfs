SET(GIDPOST_FOUND 0)

IF(EXISTS "${CFS_SOURCE_DIR}/source/DataInOut/SimInOut/GiD/gidpost/gidpost.h")
  SET(GIDPOST_FOUND 1)
ENDIF(EXISTS "${CFS_SOURCE_DIR}/source/DataInOut/SimInOut/GiD/gidpost/gidpost.h")

IF(GIDPOST_FOUND)
  #-----------------------------------------------------------------------------
  # Determine version of GiDpost by reading its version header
  #-----------------------------------------------------------------------------
  FILE(READ
    "${CFS_SOURCE_DIR}/source/DataInOut/SimInOut/GiD/gidpost/gidpost.h"
    GIDPOST_HEADER
    )

  STRING(REGEX MATCH " gidpost [0-9]+\\.[0-9]+"
    CFS_GIDPOST_VERSION
    "${GIDPOST_HEADER}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+"
        CFS_GIDPOST_VERSION "${CFS_GIDPOST_VERSION}")

#  MESSAGE("CFS_GIDPOST_VERSION ${CFS_GIDPOST_VERSION}")

ENDIF(GIDPOST_FOUND)

IF(NOT GIDPOST_FOUND)
  MESSAGE("Warning: GIDPOST could not be found! Please specify proper paths.")
ENDIF(NOT GIDPOST_FOUND)
