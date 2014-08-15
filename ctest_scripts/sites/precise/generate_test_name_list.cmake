# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  SET(${TEST_NAME_LIST}
    linux64_trunk_gcc_release 
    linux64_trunkOld_gcc_release
  )
ENDMACRO()

