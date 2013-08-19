MESSAGE("Update working copies...")

SET(UPDATE_SCRIPTS
  ctest_update_fespace_source.cmake
  ctest_update_fespace_testsuite.cmake
)

FOREACH(SCRIPT IN ITEMS ${UPDATE_SCRIPTS})

  # Checkout or update CFS++ FeSpace
  EXECUTE_PROCESS(
    COMMAND ${CTEST_COMMAND} -V -DHOME:PATH=${HOME} "-DSITE_DIR:PATH=${SITE_DIR}" -S "${SITE_DIR}/${SCRIPT}"
    RESULT_VARIABLE RETVAL
    )

  MESSAGE("RETVAL ${RETVAL}")

ENDFOREACH()
