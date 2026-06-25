# Minimal 2D FSI: a cantilever beam loaded by an analytic fluid traction supplied
# through preCICE by the mock partner (mock_fluid). cfs is the "Solid"
# participant; the coupling is serial-explicit, hence deterministic. See
# flapFSI2D.xml / make_mesh.py for the (traceable) physics.
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DEPSILON:STRING=1e-5
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DPRECICE_CONFIG:STRING=precice-config.xml
  -DMOCK_SPEC:STRING=mock_fluid.spec
  -P ${CFS_COUPLED_TEST}
)

# A coupled test runs two participants (mock_fluid + cfs). Tell ctest it uses
# 2 processes, serialise coupled tests so their preCICE m2n sockets do not clash
# under `ctest -j`, and label them so they can be selected/excluded
# (ctest -L preCICE / -LE preCICE).
set_tests_properties(${TEST_NAME} PROPERTIES
  PROCESSORS 2
  RESOURCE_LOCK precice_socket
  LABELS "preCICE")
