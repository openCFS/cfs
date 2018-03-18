# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  SET(${TEST_NAME_LIST}
    update_cfs_trunk
    update_testsuite_trunk
    linux64_trunk_gcc-6_release
    linux64_trunk_gcc-6_debug
    linux64_trunk_gcc-7_release
    linux64_trunk_gcc_release
    linux64_trunk_gcc_debug
    linux64_trunk_no-opt_gcc_release
    linux64_trunk_gcc-7_MINIMAL_release
  )
ENDMACRO()
