# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  IF(${SITE_DIR} MATCHES "master_stingl")
    SET(${TEST_NAME_LIST}
      linux64_trunk_gcc_release
      linux64_trunk_gcc_debug
    )
  ELSE()
    SET(${TEST_NAME_LIST}
      linux64_shared_opt_gcc_release
      linux64_shared_opt_gcc_mpi_release
      linux64_shared_opt_gcc_debug
      linux64_shared_opt_gcc_mpi_debug
      linux64_shared_opt_gcc_doc_doxygen
    )
  ENDIF()
ENDMACRO()

