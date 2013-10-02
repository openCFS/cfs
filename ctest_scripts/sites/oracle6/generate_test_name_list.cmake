# Manually set list of test names for the current site.
MACRO(GENERATE_TEST_NAME_LIST TEST_NAME_LIST)
  SET(${TEST_NAME_LIST}
    linux32_fespace_macosx64_release
    linux32_fespace_macosx32_release
    linux64_fespace_gcc_release
    linux64_fespace_mingw64_release
    linux64_trunk_gcc_release
    linux64_fespace_open64_release
    linux64_fespace_clang_release
  )
ENDMACRO()

