# Attention: Special flag for muParser, as gcc > 4.2 seems to produce faulty
# code with -O3, which could be tracked down to the following compiler flag
#add_definitions(-fno-inline-functions)

# prevent strict aliasing for due to bug in boost signal library v1.35
IF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
set_source_files_properties(mathParser.cc
                            PROPERTIES COMPILE_FLAGS
                            -fno-strict-aliasing)
ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")


SET(MATHPARSER_SRCS mathParser.cc registerfunc.cc)

IF(USE_OPENMP)
  SET(MATHPARSER_SRCS ${MATHPARSER_SRCS} mathParserOMP.cc)
ENDIF(USE_OPENMP)

ADD_LIBRARY(mathparser STATIC ${MATHPARSER_SRCS})

SET(TARGET_LL
  datainout
  utils
  ${MUPARSER_LIBRARY}
)

TARGET_LINK_LIBRARIES(mathparser ${TARGET_LL})

ADD_SUBDIRECTORY(pymuparser)