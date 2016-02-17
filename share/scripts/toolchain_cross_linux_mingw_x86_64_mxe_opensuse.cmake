# the name of the target operating system
SET(CMAKE_SYSTEM_NAME Windows)

SET(CMAKE_TOOLCHAIN_PREFIX "x86_64-w64-mingw32")
SET(BUILD_PREFIX "x86_64-unknown-linux-gnu")

# which compilers to use for C and C++
SET(CMAKE_C_COMPILER /opt/mxe/usr/bin/x86_64-w64-mingw32.static-gcc)
SET(CMAKE_CXX_COMPILER /opt/mxe/usr/bin/x86_64-w64-mingw32.static-g++)
SET(CMAKE_Fortran_COMPILER /opt/mxe/usr/bin/x86_64-w64-mingw32.static-gfortran)
SET(CMAKE_RC_COMPILER /opt/mxe/usr/bin/x86_64-w64-mingw32.static-windres)
SET(CMAKE_RANLIB /opt/mxe/usr/bin/x86_64-w64-mingw32.static-ranlib)

# here is the target environment located
SET(CMAKE_FIND_ROOT_PATH /opt/mxe/usr/ )

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


FIND_PROGRAM(WINE_EXECUTABLE
  NAMES wine
  )
set(TARGET_SYSTEM_EMULATOR ${WINE_EXECUTABLE})
