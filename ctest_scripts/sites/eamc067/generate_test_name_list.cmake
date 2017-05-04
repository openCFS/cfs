# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  SET(${TEST_NAME_LIST}
    linux64_mpolzer_mingw_debug
    linux64_mpolzer_gcc_debug
    linux64_mpolzer_gcc_release
    linux64_mpolzer_mingw_release
  )
ENDMACRO()

