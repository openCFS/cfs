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
#  Make a local copy of CFS++ Modelling Manual
# ===========================================================================
SET(CFS_MODELLING_MANUAL_DIR
  "${HOME}/Documents/dev/MODELLING_MANUAL_NIGHTLY")

IF(EXISTS "${CFS_MODELLING_MANUAL_DIR}")
  FILE(REMOVE_RECURSE "${CFS_MODELLING_MANUAL_DIR}")
ENDIF()

FILE(MAKE_DIRECTORY "${CFS_MODELLING_MANUAL_DIR}")

EXECUTE_PROCESS(
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${HOME}/Documents/dev/NIGHTLY/MODELLING_MANUAL_NIGHTLY" "${CFS_MODELLING_MANUAL_DIR}")

# Downgrade Graphviz to 2.28 with PNG support
# EXECUTE_PROCESS(
#  COMMAND sudo yum -y remove graphviz
# )

# EXECUTE_PROCESS(
#   COMMAND sudo yum -y install http://www.graphviz.org/pub/graphviz/stable/redhat/el6/x86_64/os/graphviz-2.28.0-1.el6.x86_64.rpm
# )

