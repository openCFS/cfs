MESSAGE("Jipi in site_specific_init.cmake on wiki!")

MESSAGE("Update working copies...")

# Checkout or update CFS++ FeSpace
EXECUTE_PROCESS(
  COMMAND ${CTEST_COMMAND} -V -DSITE_DIR:PATH=${SITE_DIR} -S "${SITE_DIR}/ctest_update_fespace_wien.cmake"
#  OUTPUT_VARIABLE DAYOFWEEK
  RESULT_VARIABLE RETVAL
  )


# Checkout or update test suite
#ctest -V -S ctest_update_fespace_testsuite_wien.cmake

# Checkout or update CFS++
#ctest -V -S ctest_update_trunk_wien.cmake

# Checkout or update test suite
#ctest -V -S ctest_update_trunk_testsuite_wien.cmake

#  Checkout or update CFSDEPS (Trunk)
#ctest -V -S ctest_update_cfsdeps_trunk_wien.cmake

