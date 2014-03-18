# http://www.cmake.org/Wiki/CMake_Cross_Compiling

# the name of the target operating system
SET(CMAKE_SYSTEM_NAME Windows)
#include(CMakeForceCompiler)
# SET(CMAKE_FIND_ROOT_PATH  /usr)

#IF("${GNU_HOST}" STREQUAL "")
#    SET(GNU_HOST i586-mingw32msvc)
#ENDIF()
# Prefix detection only works with compiler id "GNU"
#CMAKE_FORCE_C_COMPILER(x86_64-w64-mingw32-gcc GNU)
#CMAKE_FORCE_CXX_COMPILER(/usr/bin/x86_64-w64-mingw32-g++ GNU)

SET(CMAKE_TOOLCHAIN_PREFIX "x86_64-w64-mingw32")

# which compilers to use for C and C++
SET(CMAKE_C_COMPILER /usr/bin/x86_64-w64-mingw32-gcc)
SET(CMAKE_CXX_COMPILER /usr/bin/x86_64-w64-mingw32-g++)
SET(CMAKE_Fortran_COMPILER /usr/bin/x86_64-w64-mingw32-gfortran)

# CMake doesn't automatically look for prefixed 'windres', do it manually:
SET(CMAKE_RC_COMPILER /usr/bin/x86_64-w64-mingw32-windres)

SET(CMAKE_FIND_ROOT_PATH
  /usr/x86_64-w64-mingw32
  /usr/lib/gcc/x86_64-w64-mingw32
)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
