#-------------------------------------------------------------------------------
# Try to find Python
#-------------------------------------------------------------------------------
SET(PYTHON_FOUND 0)

#-------------------------------------------------------------------------------
# Determine paths of PYTHON libraries.
#-------------------------------------------------------------------------------
SET (PYTHON_POSSIBLE_LIB_PATHS
  /usr/lib64
  /usr/lib
  /usr/local/lib64
  /usr/local/lib
)

SET(PYTHON_NAMES "libpython2.7.so;libpython2.6.so;libpython2.5.a;libpython2.5.so;libpython2.4.a;libpython2.4.so;libpython2.5.dylib")

IF(NOT CFS_DISTRO STREQUAL "UBUNTU")
  SET(PYTHON_NAMES "libpython2.6.a;${PYTHON_NAMES}")
ENDIF(NOT CFS_DISTRO STREQUAL "UBUNTU")


FIND_FILE(PYTHON_LIBRARY
  NAMES ${PYTHON_NAMES}
  PATHS ${PYTHON_POSSIBLE_LIB_PATHS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  )

#-------------------------------------------------------------------------------
# Mark paths of PYTHON libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(PYTHON_LIBRARY)

#-------------------------------------------------------------------------------
# Look for PYTHON header.
#-------------------------------------------------------------------------------
SET (PYTHON_POSSIBLE_INCLUDE_PATHS
  /usr/include/python2.7
  /usr/include/python2.6
  /usr/include/python2.5
  /usr/include/python2.4
  )

FIND_PATH(PYTHON_INCLUDE_PATH
  NAMES Python.h
  PATHS ${PYTHON_POSSIBLE_INCLUDE_PATHS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  )

MARK_AS_ADVANCED(PYTHON_INCLUDE_PATH)


IF(PYTHON_LIBRARY AND PYTHON_INCLUDE_PATH)
  SET(PYTHON_FOUND 1)
ELSE(PYTHON_LIBRARY AND PYTHON_INCLUDE_PATH)
  FIND_PACKAGE(Python)
ENDIF(PYTHON_LIBRARY AND PYTHON_INCLUDE_PATH)


MARK_AS_ADVANCED(
  PYTHON_INCLUDE_PATH
  PYTHON_LIBRARY
  PYTHON_LIBRARY_DEBUG
  )

IF(PYTHON_INCLUDE_PATH)

  #-----------------------------------------------------------------------------
  # Determine version of PYTHON by preprocessing its version header
  #-----------------------------------------------------------------------------
  SET(PYTHON_TEST_FILE "${CFS_BINARY_DIR}/include/get_Python_infos.hh")
  FILE(WRITE ${PYTHON_TEST_FILE}
    "#include <${PYTHON_INCLUDE_PATH}/Python.h>\n\n"
    "<CFS_PYTHON_VERSION>PY_VERSION</CFS_PYTHON_VERSION>\n")

  IF(CMAKE_VC_COMPILER_TESTS_RUN)
    SET(PYTHON_VERSION_CMD "${CFS_CXX_COMPILER} /nologo /EP ${PYTHON_TEST_FILE}")
  ELSE(CMAKE_VC_COMPILER_TESTS_RUN)
    SET(PYTHON_VERSION_CMD "${CMAKE_CXX_COMPILER} -E ${PYTHON_TEST_FILE}")
  ENDIF(CMAKE_VC_COMPILER_TESTS_RUN)

  EXEC_PROGRAM(${PYTHON_VERSION_CMD}
    ARGS
    OUTPUT_VARIABLE PYTHON_VERSION_STR
    RETURN_VALUE RETVAL)

  FILE(REMOVE ${PYTHON_TEST_FILE})

  STRING(REGEX MATCH "<CFS_PYTHON_VERSION>.*</CFS_PYTHON_VERSION>" CFS_PYTHON_VERSION "${PYTHON_VERSION_STR}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+" CFS_PYTHON_VERSION "${CFS_PYTHON_VERSION}")

  # MESSAGE("CFS_PYTHON_VERSION ${CFS_PYTHON_VERSION}")

ENDIF(PYTHON_INCLUDE_PATH)

IF(NOT PYTHON_FOUND)
  MESSAGE("Warning: Python could not be found! Please specify proper paths.")
ENDIF(NOT PYTHON_FOUND)
