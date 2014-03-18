# ===========================================================================
#  Make a local copy of trunk test suite, since tests from different VBoxes
#  might get race conditions otherwise.
# ===========================================================================
SET(CFS_TESTSUITE_TRUNK_DIR
  "${HOME}/Documents/dev/CFS_TESTSUITE_TRUNK")

IF(EXISTS "${CFS_TESTSUITE_TRUNK_DIR}")
  FILE(REMOVE_RECURSE "${CFS_TESTSUITE_TRUNK_DIR}")
ENDIF()

FILE(MAKE_DIRECTORY "${CFS_TESTSUITE_TRUNK_DIR}")

EXECUTE_PROCESS(
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${HOME}/Documents/dev/NIGHTLY/CFS_TESTSUITE_TRUNK" "${CFS_TESTSUITE_TRUNK_DIR}")

# ===========================================================================
#  Create a proper symlink to libXmu.so which is needed by CGNS lib.
# ===========================================================================

EXECUTE_PROCESS(
  COMMAND sudo ${CMAKE_COMMAND} -E remove "/usr/lib/libXmu.so"
)

EXECUTE_PROCESS(
  COMMAND sudo ${CMAKE_COMMAND} -E create_symlink "/usr/lib/x86_64-linux-gnu/libXmu.so.6.2.0" "/usr/lib/libXmu.so"
)

