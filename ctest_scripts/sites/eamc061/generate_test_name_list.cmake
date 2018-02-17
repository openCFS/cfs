# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  IF(SITE_DIR MATCHES "trunk")
    SET(${TEST_NAME_LIST} 
     linux64_trunk_gcc_release
     linux64_trunk_gcc_debug
     linux64_trunk_gcc_serial_debug
     linux64_trunk_intel_release
     linux64_trunk_intel_debug
     linux64_trunk_clang_release
     linux64_trunk_gcc_all_release)
  ELSE()
    SET(${TEST_NAME_LIST}
     linux64_shared_opt_gcc_release
     linux64_shared_opt_clang_release
     linux64_shared_opt_intel_release
     linux64_shared_opt_gcc_debug
     linux64_shared_opt_no_prec_gcc_release
     linux64_shared_opt_mingw_release)
  ENDIF()  
ENDMACRO()
