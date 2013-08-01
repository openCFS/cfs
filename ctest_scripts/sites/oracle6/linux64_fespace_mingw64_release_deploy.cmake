#-----------------------------------------------------------------------------
# Set binary directory
#-----------------------------------------------------------------------------
SET(CTEST_BINARY_DIRECTORY "$ENV{HOME}/Documents/dev/CFS_BUILD_NIGHTLY")

message("  Searching for cfs and cfstool executables...")
SET(CFS_BIN_DIR "${CTEST_BINARY_DIRECTORY}/bin/MINGW_X86_64")
SET(CFS_BINARY "${CFS_BIN_DIR}/cfsbin.exe")
SET(CFSTOOL_BINARY "${CFS_BIN_DIR}/cfstoolbin.exe")

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
  )

  # Copy some extra .dlls from MinGW runtime
  SET(DLLS
    libgcc_s_sjlj-1.dll
    libstdc++-6.dll
    libgfortran-3.dll
    libgomp-1.dll
    pthreadGC2.dll
    pthreadGCE2.dll
  )

  FOREACH(DLL IN ITEMS ${DLLS})
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E copy "/usr/x86_64-w64-mingw32/sys-root/mingw/bin/${DLL}" "${CFS_BIN_DIR}"
      RESULT_VARIABLE RETVAL
    )
  ENDFOREACH()

  message("  Packing binaries...")
  EXECUTE_PROCESS(
    COMMAND zip -yr "${NIGHTLY_ARCHIVES_DIR}/${HOSTNAME}_${TEST_NAME}.zip" CFS_BUILD_NIGHTLY
    WORKING_DIRECTORY "${CTEST_BINARY_DIRECTORY}/.."
    RESULT_VARIABLE RETVAL
  )

ENDIF()
