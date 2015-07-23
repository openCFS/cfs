#-----------------------------------------------------------------------------
# Set binary directory
#-----------------------------------------------------------------------------
SET(CTEST_BINARY_DIRECTORY "$ENV{HOME}/Documents/dev/CFS_BUILD_NIGHTLY")

message("  Searching for cfs and cfstool executables...")
SET(CFS_BIN_DIR "${CTEST_BINARY_DIRECTORY}/bin/MINGW_X86_64")
SET(CFS_BINARY "${CFS_BIN_DIR}/cfsbin")
SET(CFSTOOL_BINARY "${CFS_BIN_DIR}/cfstoolbin")

MESSAGE("CFS_BINARY ${CFS_BINARY}")
MESSAGE("CFSTOOL_BINARY ${CFSTOOL_BINARY}")

IF(EXISTS "${CFS_BINARY}" AND
   EXISTS "${CFSTOOL_BINARY}")

  message("  Cleaning up binary directory...")

  # Remove unnecessary subdirs before packing archive.
  FILE(REMOVE_RECURSE
    "${CTEST_BINARY_DIRECTORY}/cfsdeps"
    "${CTEST_BINARY_DIRECTORY}/include"
    "${CTEST_BINARY_DIRECTORY}/tmp"
    "${CTEST_BINARY_DIRECTORY}/testsuite"
    "${CTEST_BINARY_DIRECTORY}/source"
    "${CTEST_BINARY_DIRECTORY}/CMakeFiles"
    "${CTEST_BINARY_DIRECTORY}/Testing"
  )

  file(GLOB_RECURSE STATIC_ARCHIVES "${CTEST_BINARY_DIRECTORY}/lib64/*.a")
  FILE(REMOVE ${STATIC_ARCHIVES})

  message("  Packing binaries...")
  EXECUTE_PROCESS(
    COMMAND zip -yr "${NIGHTLY_ARCHIVES_DIR}/${HOSTNAME}_${TEST_NAME}.zip" CFS_BUILD_NIGHTLY
    WORKING_DIRECTORY "${CTEST_BINARY_DIRECTORY}/.."
    RESULT_VARIABLE RETVAL
  )

ENDIF()
