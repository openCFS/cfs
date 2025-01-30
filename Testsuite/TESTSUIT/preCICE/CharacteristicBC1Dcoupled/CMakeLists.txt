# openCFS<->openCFS impedance-matched (characteristic) acoustic coupling. A plane-wave
# duct is split at x=Lf into two openCFS acoustic participants coupled via preCICE: the
# "Fluid" partner (fluid.xml) launches a 1-cycle 5 kHz burst and exchanges the characteristic
# invariants at the interface with the "Acoustic" main participant (this test, compared
# against the .h5ref). Replaces the original OpenFOAM Fluid domain with a second openCFS
# instance, so it needs the cfs+cfs coupled-run driver (CoupledCfsTest.cmake). Serial-explicit
# => deterministic. See CharacteristicBC1Dcoupled.xml / fluid.xml / precice-config.xml / make_mesh.py.
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DEPSILON:STRING=1e-4
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DPARTNER_BASENAME:STRING=fluid
  -DPRECICE_CONFIG:STRING=precice-config.xml
  # Absolute difference: the acoustic domain is quiescent until the burst crosses
  # S_couple from the Fluid participant, so the first steps are machine noise
  # (< 1e-11) and a relative comparison of them is solver/platform noise. 1e-4
  # absolute still locks the physical transmission tightly.
  -DCFSTOOL_MODE=absL2diff
  -P ${CFS_COUPLED_CFS_TEST}
)

# A coupled test runs two participants (two cfs instances). Tell ctest it uses
# 2 processes, serialise coupled tests so their preCICE m2n sockets do not clash
# under `ctest -j`, and label them so they can be selected/excluded
# (ctest -L preCICE / -LE preCICE).
set_tests_properties(${TEST_NAME} PROPERTIES
  PROCESSORS 2
  RESOURCE_LOCK precice_socket
  LABELS "preCICE")
