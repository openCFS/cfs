#-------------------------------------------------------------------------------
# Try to find Tcl with CMake standard package.
#-------------------------------------------------------------------------------
SET(TCL_FOUND 0)

FIND_PACKAGE(TCL)

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
