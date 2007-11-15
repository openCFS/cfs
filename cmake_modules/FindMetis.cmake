SET(METIS_FOUND 0)

#-------------------------------------------------------------------------------
# Determine path of METIS library.
#-------------------------------------------------------------------------------
FIND_LIBRARY(METIS_LIBRARY
  NAMES metis
  PATHS ${CFSDEPS_LIBRARY_DIR}
  )

#-------------------------------------------------------------------------------
# Mark path of METIS library as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(METIS_LIBRARY)

FIND_PATH(METIS_INCLUDE_DIR metis.h 
  ${CFSDEPS_INCLUDE_DIR} 
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH)

MARK_AS_ADVANCED(METIS_INCLUDE_DIR)


IF(METIS_LIBRARY AND METIS_INCLUDE_DIR)
  SET(METIS_FOUND 1)
ENDIF(METIS_LIBRARY AND METIS_INCLUDE_DIR)

IF(METIS_FOUND)

  #-----------------------------------------------------------------------------
  # Determine version of METIS by preprocessing its version header
  #-----------------------------------------------------------------------------
  SET(METIS_TEST_FILE "${CFS_BINARY_DIR}/include/get_metis_infos.hh")
  FILE(WRITE ${METIS_TEST_FILE}
    "#include <${CFSDEPS_INCLUDE_DIR}/defs.h>\n\n"
    "<CFS_METIS_VERSION>METISTITLE</CFS_METIS_VERSION>\n")

  SET(METIS_VERSION_CMD "${CMAKE_CXX_COMPILER} -E ${METIS_TEST_FILE}")

  EXEC_PROGRAM(${METIS_VERSION_CMD}
    ARGS
    OUTPUT_VARIABLE METIS_VERSION_STR
    RETURN_VALUE RETVAL)

  #  MESSAGE("METIS_VERSION_STR ${METIS_VERSION_STR}")
  FILE(REMOVE ${METIS_TEST_FILE})

  STRING(REGEX MATCH
    "<CFS_METIS_VERSION>.*</CFS_METIS_VERSION>"
    CFS_METIS_VERSION "${METIS_VERSION_STR}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
    CFS_METIS_VERSION "${CFS_METIS_VERSION}")

#  MESSAGE("CFS_METIS_VERSION ${CFS_METIS_VERSION}")

#  MARK_AS_ADVANCED(METIS_LIBRARY)

ENDIF(METIS_FOUND)

IF(NOT METIS_FOUND)
  MESSAGE("Warning: METIS could not be found! Please specify proper paths.")
ENDIF(NOT METIS_FOUND)
