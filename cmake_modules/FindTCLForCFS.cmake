#-------------------------------------------------------------------------------
# Try to find Tcl with CMake standard package.
#-------------------------------------------------------------------------------
SET(TCL_FOUND 0)

#-------------------------------------------------------------------------------
# Determine paths of TCL libraries.
#-------------------------------------------------------------------------------
SET (TCL_POSSIBLE_LIB_PATHS
  ${CFSDEPS_LIBRARY_DIR}
  /usr/lib64
  /usr/lib
  /usr/local/lib64
  /usr/local/lib
)

FIND_LIBRARY(TCL_LIBRARY
  NAMES tcl8.4
  PATHS ${TCL_POSSIBLE_LIB_PATHS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  )

#-------------------------------------------------------------------------------
# Mark paths of TCL libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(TCL_LIBRARY)


#-------------------------------------------------------------------------------
# Look for TCL header.
#-------------------------------------------------------------------------------
SET (TCL_POSSIBLE_INCLUDE_PATHS
  ${CFSDEPS_INCLUDE_DIR}
  /usr/include
  /usr/local/include
  )

FIND_PATH(TCL_INCLUDE_PATH
  NAMES tcl.h
  PATHS ${TCL_POSSIBLE_INCLUDE_PATHS}
  )

MARK_AS_ADVANCED(TCL_INCLUDE_PATH)


IF(TCL_LIBRARY AND XERCES_INCLUDE_PATH)
  SET(TCL_FOUND 1)
ELSE(TCL_LIBRARY AND XERCES_INCLUDE_PATH)
  FIND_PACKAGE(TCL)
ENDIF(TCL_LIBRARY AND XERCES_INCLUDE_PATH)


MARK_AS_ADVANCED(
  TCL_INCLUDE_PATH
  TCL_LIBRARY
  TCL_LIBRARY_DEBUG
  TCL_TCLSH
  TK_INCLUDE_PATH
  TK_LIBRARY
  TK_LIBRARY_DEBUG
  TK_WISH
  )

IF(TCL_INCLUDE_PATH)

  #-----------------------------------------------------------------------------
  # Determine version of TCL by preprocessing its version header
  #-----------------------------------------------------------------------------
  SET(TCL_TEST_FILE "${CFS_BINARY_DIR}/include/get_tcl_infos.hh")
  FILE(WRITE ${TCL_TEST_FILE}
    "#include <${TCL_INCLUDE_PATH}/tcl.h>\n\n"
    "<CFS_TCL_VERSION>TCL_PATCH_LEVEL</CFS_TCL_VERSION>\n")

  IF(CMAKE_VC_COMPILER_TESTS_RUN)
    SET(TCL_VERSION_CMD "${CFS_CXX_COMPILER} /nologo /EP ${TCL_TEST_FILE}")
  ELSE(CMAKE_VC_COMPILER_TESTS_RUN)
    SET(TCL_VERSION_CMD "${CMAKE_CXX_COMPILER} -E ${TCL_TEST_FILE}")
  ENDIF(CMAKE_VC_COMPILER_TESTS_RUN)

  EXEC_PROGRAM(${TCL_VERSION_CMD}
    ARGS
    OUTPUT_VARIABLE TCL_VERSION_STR
    RETURN_VALUE RETVAL)

  FILE(REMOVE ${TCL_TEST_FILE})

  STRING(REGEX MATCH "<CFS_TCL_VERSION>.*</CFS_TCL_VERSION>" CFS_TCL_VERSION "${TCL_VERSION_STR}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" CFS_TCL_VERSION "${CFS_TCL_VERSION}")

#  MESSAGE("CFS_TCL_VERSION ${CFS_TCL_VERSION}")

  SET(USE_SCRIPTING 1)
  SET(USE_SCRIPTING_TCL 1)
  SET(TCL_FOUND 1)

ENDIF(TCL_INCLUDE_PATH)

IF(NOT TCL_FOUND)
  MESSAGE("Warning: TCL could not be found! Please specify proper paths.")
ENDIF(NOT TCL_FOUND)
