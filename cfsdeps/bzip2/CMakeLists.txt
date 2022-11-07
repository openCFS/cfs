CMAKE_MINIMUM_REQUIRED(VERSION 3.6)

PROJECT(bzip2 C)

set(INSTALL_LIB_DIR
  "${CMAKE_INSTALL_PREFIX}/${LIB_SUFFIX}"
  CACHE PATH "Installation directory for libraries")
set(INSTALL_INC_DIR
  "${CMAKE_INSTALL_PREFIX}/include"
  CACHE PATH "Installation directory for headers")

SET(BZIP2_SRCS
  blocksort.c
  huffman.c
  crctable.c
  randtable.c
  compress.c
  decompress.c
  bzlib.c)

SET(BZIP2_LIBNAME bz2)

ADD_DEFINITIONS(-D_FILE_OFFSET_BITS=64)

IF(WIN32) 
  IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    SET(BZIP2_LIBNAME bz2d)
    ADD_DEFINITIONS(-UBZ_UNIX -DBZ_LCCWIN32=1)
  ENDIF()

  IF(MSVC)
    ADD_DEFINITIONS(/D_CRT_SECURE_NO_WARNINGS)
  ENDIF()
ENDIF(WIN32)
IF(UNIX)
  ADD_DEFINITIONS(-DBZ_UNIX=1 -UBZ_LCCWIN32)
ENDIF(UNIX)

# MESSAGE(FATAL_ERROR "${BZIP2_LIBNAME}")
ADD_LIBRARY(${BZIP2_LIBNAME} ${BZIP2_SRCS})

# parse the full version number from zlib.h and include in ZLIB_FULL_VERSION
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/bzlib_private.h _bzlib_private_h_contents)
string(REGEX REPLACE ".*#define[ \t]+BZ_VERSION[ \t]+\"([0-9A-Za-z.]+).*"
  "\\1" BZIP2_FULL_VERSION ${_bzlib_private_h_contents})

set_target_properties(${BZIP2_LIBNAME} PROPERTIES SOVERSION 1)

# This property causes shared libraries on Linux to have the full version
# encoded into their final filename.  

# This has no effect with MSVC, on that platform the version info for
# the DLL comes from the resource file win32/zlib1.rc
set_target_properties(${BZIP2_LIBNAME} PROPERTIES VERSION ${BZIP2_FULL_VERSION})

INSTALL(TARGETS ${BZIP2_LIBNAME}
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}")

INSTALL(FILES bzlib.h DESTINATION include)