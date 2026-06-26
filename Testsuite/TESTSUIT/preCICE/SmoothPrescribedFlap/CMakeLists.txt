# SmoothPDE "prescribed displacement" mode fed via preCICE. The mock partner
# (mock_fluid) plays "Dummy": it writes a uniform, linearly-ramped DisplacementVolume
# on a single vertex; cfs is the "Solid" participant and injects the received field
# directly into SMOOTH_DISPLACEMENT (no solve), exporting smoothDisplacement. With
# nearest-neighbour mapping every node gets the same vector, so the result is the
# deterministic ramp u(t) = (1e-3, 5e-4, 0)*min(t/0.1, 1). Serial-explicit coupling
# => deterministic. See SmoothPrescribedFlap.xml / precice-config.xml / make_mesh.py.
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
