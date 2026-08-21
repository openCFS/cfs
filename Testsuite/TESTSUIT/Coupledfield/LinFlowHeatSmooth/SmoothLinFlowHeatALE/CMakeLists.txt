#-------------------------------------------------------------------------------
# SmoothPDE -> LinFlowPDE + HeatPDE ALE function test
#
# SmoothPDE generates
#
#     v_g = (0.1, 0) m/s
#
# and supplies smoothVelocity to LinFlowPDE and HeatPDE.
#
# Both PDEs use the identical movingMeshId "moveGridID".
#
# The complete transient solution was independently validated against
# an otherwise identical LinFlow-Heat simulation with directly
# prescribed grid velocity.
#
# Pressure, LinFlow velocity and temperature agree to approximately
# 1e-13 relative L2 error.
#
# SmoothPDE results show small OpenMP-dependent numerical deviations.
# Therefore the regression tolerance is set to 1e-6.
#-------------------------------------------------------------------------------


GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")


ADD_TEST(
  ${TEST_NAME}

  ${CMAKE_COMMAND}

  -DEPSILON:STRING=1e-6

  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}

  -DADD_TEST_FILES:STRING=CavityThermoviscous.h5ref

  -P ${CFS_STANDARD_TEST}
)


SET_TESTS_PROPERTIES(
  ${TEST_NAME}

  PROPERTIES

  LABELS "coupled;linflow;heat;smooth"
)