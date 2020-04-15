# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  if(SITE_DIR MATCHES "master")
    set(${TEST_NAME_LIST} 
     linux64_master_gcc_release
     linux64_master_gcc_debug
     linux64_master_clang_release
  else()
    set(${TEST_NAME_LIST}
     linux64_shared_opt_gcc_release
     linux64_shared_opt_gcc_debug
     linux64_shared_opt_no_prec_gcc_release
     linux64_shared_opt_openblas_release

     linux64_shared_opt_gcc8_release
     # linux64_shared_opt_intel_release
     linux64_shared_opt_clang_release)
  endif()  
ENDMACRO()
