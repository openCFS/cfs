#-------------------------------------------------------------------------------
# SmoothPDE -> HeatPDE ALE function test
#
# SmoothPDE generates a constant grid velocity
#
#     v_g = (-0.1, 0) m/s
#
# which is passed to HeatPDE through smoothVelocity and movingMeshList.
#
# The setup was independently validated against an otherwise identical
# HeatPDE simulation with directly prescribed grid velocity.
#
# Temperature differences between both formulations were below 1e-12
# for all 10 transient time steps.
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------

GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")


#-------------------------------------------------------------------------------
# Define the regression test.
#-------------------------------------------------------------------------------

ADD_TEST(${TEST_NAME}

  ${CMAKE_COMMAND}

  -DEPSILON:STRING=1e-12

  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}

  -P ${CFS_STANDARD_TEST}

)


#-------------------------------------------------------------------------------
# Test labels
#-------------------------------------------------------------------------------

SET_TESTS_PROPERTIES(
  ${TEST_NAME}
  PROPERTIES
  LABELS "coupled;heat;smooth"
)
