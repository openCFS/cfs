SET(CGAL_FOUND 0)

#-------------------------------------------------------------------------------
# Determine path of CGAL library.
#-------------------------------------------------------------------------------
FIND_LIBRARY(CGAL_LIBRARY
  NAMES CGAL
  PATHS ${CFSDEPS_LIBRARY_DIR}
  )

#-------------------------------------------------------------------------------
# Mark path of CGAL library as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(CGAL_LIBRARY)

FIND_PATH(CGAL_INCLUDE_DIR CGAL/version.h 
  ${CFSDEPS_INCLUDE_DIR} 
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH)

MARK_AS_ADVANCED(CGAL_INCLUDE_DIR)


IF(CGAL_LIBRARY AND CGAL_INCLUDE_DIR)
  SET(CGAL_FOUND 1)
ENDIF(CGAL_LIBRARY AND CGAL_INCLUDE_DIR)

IF(CGAL_FOUND)

  #-----------------------------------------------------------------------------
  # Determine version of CGAL by preprocessing its version header
  #-----------------------------------------------------------------------------
  SET(CGAL_TEST_FILE "${CFS_BINARY_DIR}/include/get_cgal_infos.hh")
  FILE(WRITE ${CGAL_TEST_FILE}
    "#include <${CGAL_INCLUDE_DIR}/CGAL/version.h>\n\n"
    "<CFS_CGAL_VERSION>CGAL_VERSION</CFS_CGAL_VERSION>\n")

  SET(CGAL_VERSION_CMD "${CMAKE_CXX_COMPILER} -E ${CGAL_TEST_FILE}")

  EXEC_PROGRAM(${CGAL_VERSION_CMD}
    ARGS
    OUTPUT_VARIABLE CGAL_VERSION_STR
    RETURN_VALUE RETVAL)

#  MESSAGE("CGAL_VERSION_STR ${CGAL_VERSION_STR}")
  FILE(REMOVE ${CGAL_TEST_FILE})

  STRING(REGEX MATCH
    "<CFS_CGAL_VERSION>.*</CFS_CGAL_VERSION>"
    CFS_CGAL_VERSION "${CGAL_VERSION_STR}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+"
    CFS_CGAL_VERSION "${CFS_CGAL_VERSION}")

#  MESSAGE("CFS_CGAL_VERSION ${CFS_CGAL_VERSION}")

ENDIF(CGAL_FOUND)

IF(NOT CGAL_FOUND)
  MESSAGE("Warning: CGAL could not be found! Please specify proper paths.")
ENDIF(NOT CGAL_FOUND)
