CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(bzip2 C)

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

IF(WIN32 AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  SET(BZIP2_LIBNAME bz2d)
  ADD_DEFINITIONS(-UBZ_UNIX -DBZ_LCCWIN32=1)
ENDIF(WIN32 AND CMAKE_BUILD_TYPE STREQUAL "Debug")

IF(UNIX)
  ADD_DEFINITIONS(-DBZ_UNIX=1 -UBZ_LCCWIN32)
ENDIF(UNIX)

# MESSAGE(FATAL_ERROR "${BZIP2_LIBNAME}")
ADD_LIBRARY(${BZIP2_LIBNAME} ${BZIP2_SRCS})

IF(NOT BUILD_SHARED_LIBS)
  ADD_EXECUTABLE(bzip2 bzip2.c)
  TARGET_LINK_LIBRARIES(bzip2 ${BZIP2_LIBNAME})
  
  ADD_EXECUTABLE(bzip2recover bzip2recover.c)
  TARGET_LINK_LIBRARIES(bzip2recover ${BZIP2_LIBNAME})

  INSTALL(TARGETS bzip2 bzip2recover
    RUNTIME DESTINATION "bin/${CFS_ARCH_STR}")

  INSTALL(TARGETS ${BZIP2_LIBNAME}
    ARCHIVE DESTINATION "${LIB_SUFFIX}/${CFS_ARCH_STR}")
ELSE()

  # parse the full version number from zlib.h and include in ZLIB_FULL_VERSION
  file(READ ${CMAKE_CURRENT_SOURCE_DIR}/bzlib_private.h _bzlib_private_h_contents)
  string(REGEX REPLACE ".*#define[ \t]+BZ_VERSION[ \t]+\"([0-9A-Za-z.]+).*"
    "\\1" BZIP2_FULL_VERSION ${_bzlib_private_h_contents})

  set_target_properties(${BZIP2_LIBNAME} PROPERTIES SOVERSION 1)

  if(NOT CYGWIN)
    # This property causes shared libraries on Linux to have the full version
    # encoded into their final filename.  We disable this on Cygwin because
    # it causes cygz-${ZLIB_FULL_VERSION}.dll to be created when cygz.dll
    # seems to be the default.
    #
    # This has no effect with MSVC, on that platform the version info for
    # the DLL comes from the resource file win32/zlib1.rc
    set_target_properties(${BZIP2_LIBNAME} PROPERTIES VERSION ${BZIP2_FULL_VERSION})
  endif()

  IF(WIN32)
    INSTALL(TARGETS ${BZIP2_LIBNAME}
      ARCHIVE DESTINATION "${LIB_SUFFIX}/${CFS_ARCH_STR}"
      RUNTIME DESTINATION "bin/${CFS_ARCH_STR}")
  ELSE()
    INSTALL(TARGETS ${BZIP2_LIBNAME}
      LIBRARY DESTINATION "${LIB_SUFFIX}/${CFS_ARCH_STR}")
  ENDIF()
ENDIF(NOT BUILD_SHARED_LIBS)



INSTALL(FILES
  bzlib.h
  DESTINATION include)