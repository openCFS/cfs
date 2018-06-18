# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  IF(SITE_DIR MATCHES "simon")
    SET(${TEST_NAME_LIST} 
     linux64_simon_gcc_release
     # linux64_simon_gcc_debug
     # linux64_simon_intel_release
     # linux64_simon_intel_debug
     )
  ELSE()
    SET(${TEST_NAME_LIST}
     # linux64_shared_opt_gcc_release
     # linux64_shared_opt_gcc_debug
     # linux64_shared_opt_intel_release
     # linux64_shared_opt_intel_debug
     )
  ENDIF()  
ENDMACRO()
