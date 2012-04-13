# The muparser external project

set(muparser_source  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/muparser")
set(muparser_binary  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/muparser-build")
set(muparser_install  "${CMAKE_CURRENT_BINARY_DIR}")

FILE(MAKE_DIRECTORY "${muparser_source}")

SET(PFN "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/patches/muparser_patch.cmake")
FILE(WRITE ${PFN} "SET(muparser_source  \"${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/muparser\")\n")
FILE(APPEND ${PFN} "SET(CFS_SOURCE_DIR \"${CFS_SOURCE_DIR}\")\n\n")

FILE(APPEND ${PFN} "SET(files \"include/muParserDef.h\"\n")
FILE(APPEND ${PFN} "	  \"include/muParserBase.h\"\n")
FILE(APPEND ${PFN} "	  \"include/muParserCallback.h\"\n")
FILE(APPEND ${PFN} "	  \"include/muParserTemplateMagic.h\"\n")
FILE(APPEND ${PFN} "      \"src/muParser.cpp\"\n")
FILE(APPEND ${PFN} "	  \"src/muParserBase.cpp\"\n")
FILE(APPEND ${PFN} "	  \"src/muParserCallback.cpp\")\n\n")

FILE(APPEND ${PFN} "foreach(file IN LISTS files)\n")
FILE(APPEND ${PFN} "  MESSAGE(\"UNIX newlines for file: \${file}\")\n")
FILE(APPEND ${PFN} "  configure_file(\${muparser_source}/\${file} \${muparser_source}/dummy @ONLY NEWLINESTYLE UNIX)\n")
FILE(APPEND ${PFN} "  configure_file(\${muparser_source}/dummy \${muparser_source}/\${file} @ONLY NEWLINESTYLE UNIX)\n")
FILE(APPEND ${PFN} "endforeach(file)\n\n")

FILE(APPEND ${PFN} "EXECUTE_PROCESS(\n") 
FILE(APPEND ${PFN} "  COMMAND \${CMAKE_COMMAND} -E copy_if_different \"\${CFS_SOURCE_DIR}/cfsdeps/muparser/CMakeLists.txt\" \"\${muparser_source}/CMakeLists.txt\"\n") 
FILE(APPEND ${PFN} "  COMMAND patch -p0 -i \"\${CFS_SOURCE_DIR}/cfsdeps/muparser/muparser-strfun45.patch\"\n") 
FILE(APPEND ${PFN} "  COMMAND patch -p0 -i \"\${CFS_SOURCE_DIR}/cfsdeps/muparser/muparser-dummy-math.patch\"\n") 

# The following patch is due to some missing intrinsic functions
# in Intel compiler. Cf. Intel® C++ Compiler for Linux 
# Intrinsics Reference
# http://download.intel.com/support/performancetools/c/linux/v9/intref_cls.pdf
IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
  FILE(APPEND ${PFN} "  COMMAND patch -p0 -i \"\${CFS_SOURCE_DIR}/cfsdeps/muparser/muparser-intel-atomic.patch\"\n") 
ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")

FILE(APPEND ${PFN} ")\n") 

ExternalProject_Add(muparser
  DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cache
  URL ${MUPARSER_URL}/${MUPARSER_ZIP}
  URL_MD5 ${MUPARSER_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  SOURCE_DIR ${muparser_source}
  BINARY_DIR ${muparser_binary}
  CMAKE_ARGS
#     ${muparser_lib_args}
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
#    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG:PATH=${muparser_binary}/lib
#    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=${muparser_binary}/lib
#    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE:PATH=${muparser_binary}/lib
#    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=${muparser_binary}/lib
    -DCMAKE_INSTALL_PREFIX:PATH=${muparser_install}
    -DARCH_STR:STRING=${CFS_ARCH_STR}
#    -DARCH:STRING=${CFS_ARCH}
#    -DREV:STRING=${CFS_REV}
#    -DDIST:STRING=${CFS_DIST}
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
    ${Muparser_EXTRA_ARGS}
  )

#  ExternalProject_Add_Step(Muparser InstallUUIDInclude
#    COMMAND ${CMAKE_COMMAND} -E copy_directory ${muparser_source}/src/uuid/include/muparser/uuid ${muparser_install}/include/muparser/uuid
#    DEPENDEES install
#  )

SET(MUPARSER_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libmuparser.a CACHE FILEPATH "muparser library")
SET(MUPARSER_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "muparser include directory")

MARK_AS_ADVANCED(MUPARSER_LIBRARY)
MARK_AS_ADVANCED(MUPARSER_INCLUDE_DIR)
