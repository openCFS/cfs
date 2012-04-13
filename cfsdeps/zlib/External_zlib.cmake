# The zlib external project

set(zlib_source  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/zlib")
set(zlib_binary  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/zlib-build")
set(zlib_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(PFN "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/patches/zlib_patch.cmake")
FILE(WRITE ${PFN} "SET(zlib_source  \"${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/zlib\")\n")
FILE(APPEND ${PFN} "SET(CFS_SOURCE_DIR \"${CFS_SOURCE_DIR}\")\n\n") 
FILE(APPEND ${PFN} "EXECUTE_PROCESS(\n") 
FILE(APPEND ${PFN} "  COMMAND \${CMAKE_COMMAND} -E remove -f \"\${zlib_source}/zconf.h\"\n") 
FILE(APPEND ${PFN} "  COMMAND patch -p0 -i \"\${CFS_SOURCE_DIR}/cfsdeps/zlib/zlib-cmakelists.patch\"\n") 
FILE(APPEND ${PFN} ")\n") 

ExternalProject_Add(zlib-shared
  DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cache
  URL ${ZLIB_URL}/${ZLIB_GZ}
  URL_MD5 ${ZLIB_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  SOURCE_DIR ${zlib_source}
  BINARY_DIR ${zlib_binary}
  CMAKE_ARGS
#     ${zlib_lib_args}
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
#    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG:PATH=${zlib_binary}/lib
#    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=${zlib_binary}/lib
#    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE:PATH=${zlib_binary}/lib
#    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=${zlib_binary}/lib
    -DCMAKE_INSTALL_PREFIX:PATH=${zlib_install}
    -DARCH_STR:STRING=${CFS_ARCH_STR}
#    -DARCH:STRING=${CFS_ARCH}
#    -DREV:STRING=${CFS_REV}
#    -DDIST:STRING=${CFS_DIST}
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
    ${Zlib_EXTRA_ARGS}
  )


ExternalProject_Add(zlib-static
  DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cache
  URL ${ZLIB_URL}/${ZLIB_GZ}
  URL_MD5 ${ZLIB_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  SOURCE_DIR ${zlib_source}
  BINARY_DIR ${zlib_binary}
  CMAKE_ARGS
#     ${zlib_lib_args}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
#    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG:PATH=${zlib_binary}/lib
#    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=${zlib_binary}/lib
#    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE:PATH=${zlib_binary}/lib
#    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=${zlib_binary}/lib
    -DCMAKE_INSTALL_PREFIX:PATH=${zlib_install}
    -DARCH_STR:STRING=${CFS_ARCH_STR}
#    -DARCH:STRING=${CFS_ARCH}
#    -DREV:STRING=${CFS_REV}
#    -DDIST:STRING=${CFS_DIST}
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
    ${Zlib_EXTRA_ARGS}
  )

#  ExternalProject_Add_Step(Zlib InstallUUIDInclude
#    COMMAND ${CMAKE_COMMAND} -E copy_directory ${zlib_source}/src/uuid/include/zlib/uuid ${zlib_install}/include/zlib/uuid
#    DEPENDEES install
#  )

SET(ZLIB_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libz.a CACHE FILEPATH "zlib library")
SET(ZLIB_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "zlib include directory")

MARK_AS_ADVANCED(ZLIB_LIBRARY)
MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR)
