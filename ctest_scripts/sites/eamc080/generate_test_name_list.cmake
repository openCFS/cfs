# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  IF(SITE_DIR MATCHES "trunk")
    SET(${TEST_NAME_LIST}
      linux64_trunk_gcc_debug.ctest
      linux64_trunk_gcc_release.ctest
    )
  ELSE()
    SET(${TEST_NAME_LIST}
      linux64_shared_opt_gcc_release
      linux64_shared_opt_gcc_debug
      linux64_shared_opt_gcc_doc_doxygen
    )
  ENDIF()
ENDMACRO()

