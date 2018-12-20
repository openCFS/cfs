# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  if(SITE_DIR MATCHES "master_wien")
    set(${TEST_NAME_LIST} 
     linux64_trunk_gcc_release
     #linux64_trunk_gcc_all_release
     linux64_trunk_gcc_debug
     #linux64_trunk_gcc_serial_debug
     linux64_trunk_intel_release)
     #linux64_trunk_intel_debug
     #linux64_trunk_clang_release
     #linux64_trunk_gcc_all_release)
  elseif(SITE_DIR MATCHES "master_stingl")     
    set(${TEST_NAME_LIST} 
     linux64_trunk_gcc_release
     linux64_trunk_gcc_debug
     linux64_trunk_intel_release)
  else()
    SET(${TEST_NAME_LIST}
     linux64_shared_opt_gcc_release
     linux64_shared_opt_clang_release
     linux64_shared_opt_intel_release
     linux64_shared_opt_gcc_debug
     linux64_shared_opt_gcc_serial_debug  
     linux64_shared_opt_no_prec_gcc_release
     #linux64_shared_opt_mingw_release
     linux64_shared_opt_openblas_release)
  endif()  
ENDMACRO()
