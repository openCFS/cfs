# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  SET(${TEST_NAME_LIST}
#    linux32_trunk_macosx64_release
#    linux32_trunk_macosx32_release
    linux64_trunk_gcc_release
#    linux64_trunk_mingw64_release
#    linux64_trunkOld_gcc_release
#    linux64_trunk_open64_release
#    linux64_trunk_clang_release
  )
ENDMACRO()

