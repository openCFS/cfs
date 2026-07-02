# Full coupling combination (OpenFOAM replaced by the mock_fluid partner): iterative
# (parallel-implicit + IQN-ILS) FSI plus a one-way SmoothPDE volume deformation, with TWO
# PDEs on the cfs "Solid" participant. The mock plays "Fluid": it writes the interface Stress
# and a whole-volume DisplacementVolume and reads back the displacement; the implicit scheme
# re-iterates each window (exercising the checkpoint path in the mock + adapter). cfs reads the
# traction (mechNormalStressPrecice) into the mechanic PDE and the volume deformation
# (smoothDisplacement) into the prescribed SmoothPDE. Serial determinism: the mock writes
# analytic (time-only) data and restores its state on re-iteration. See PerpendicularFlapSmooth.xml
# / precice-config.xml / make_mesh.py.
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
