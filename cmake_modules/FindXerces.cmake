# - Find Xerces-C library.
# This module finds if Xerces is installed and determines where the
# library is. It also determines what the name of the library is.
# This code sets the following variables:
#  XERCES_LIBRARY        = path to XERCES library (xerces-c)
#  XERCES_INCLUDE_DIR    = path to XERCES headers
#  XERCES_FOUND          = bool which says if BLAS has been found
#  XERCES_VERSION_MAJOR
#  XERCES_VERSION_MINOR
#  XERCES_VERSION_REVISION
#  XERCES_VERSION

SET (XERCES_POSSIBLE_LIB_PATHS
  /usr/lib64
  /usr/lib
  /usr/local/lib64
  /usr/local/lib
)

FIND_LIBRARY(XERCES_LIBRARY
  NAMES xerces-c
  PATHS ${XERCES_POSSIBLE_LIB_PATHS}
)

SET (XERCES_POSSIBLE_INCLUDE_PATHS
  /usr/include
  /usr/local/include
)

FIND_PATH(XERCES_INCLUDE_DIR xercesc/util/XercesVersion.hpp 
  ${XERCES_POSSIBLE_INCLUDE_PATHS} )

IF(XERCES_INCLUDE_DIR)
  EXEC_PROGRAM("cat ${XERCES_INCLUDE_DIR}/xercesc/util/XercesVersion.hpp | grep '#define XERCES_VERSION_MAJOR' | tail -1 | cut -d' ' -f 3"
    ARGS
    OUTPUT_VARIABLE XERCES_VERSION_MAJOR
    RETURN_VALUE RETVAL)

  EXEC_PROGRAM("cat ${XERCES_INCLUDE_DIR}/xercesc/util/XercesVersion.hpp | grep '#define XERCES_VERSION_MINOR' | tail -1 | cut -d' ' -f 3"
    ARGS
    OUTPUT_VARIABLE XERCES_VERSION_MINOR
    RETURN_VALUE RETVAL)

  EXEC_PROGRAM("cat ${XERCES_INCLUDE_DIR}/xercesc/util/XercesVersion.hpp | grep '#define XERCES_VERSION_REVISION' | tail -1 | cut -d' ' -f 3"
    ARGS
    OUTPUT_VARIABLE XERCES_VERSION_REVISION
    RETURN_VALUE RETVAL)

  SET(XERCES_VERSION "${XERCES_VERSION_MAJOR}.${XERCES_VERSION_MINOR}.${XERCES_VERSION_REVISION}")
ENDIF(XERCES_INCLUDE_DIR)

IF(XERCES_LIBRARY)
  IF(XERCES_INCLUDE_DIR)
    SET(XERCES_FOUND 1)
#    MESSAGE("XERCES_VERSION ${XERCES_VERSION}")
  ENDIF(XERCES_INCLUDE_DIR)
ENDIF(XERCES_LIBRARY)

MARK_AS_ADVANCED(
  XERCES_LIBRARY
  XERCES_INCLUDE_DIR
  XERCES_FOUND
  XERCES_VERSION_MAJOR
  XERCES_VERSION_MINOR
  XERCES_VERSION_REVISION
  XERCES_VERSION
  )
