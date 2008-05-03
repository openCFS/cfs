# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Attention: Special flag for muParser, as gcc > 4.2 seems to produce faulty code with -O3, which could be tracked down to the following compiler flag
add_definitions(-fno-inline-functions)

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()
SET(MATHPARSER_SRCS mathParser.cc muParserBase.cpp muParser.cpp
  muParserBytecode.cpp muParserCallback.cpp muParserError.cpp
  muParserTokenReader.cpp muParserInt.cpp)

ADD_LIBRARY(mathparser STATIC ${MATHPARSER_SRCS})

SET(TARGET_LL
  utils
)

TARGET_LINK_LIBRARIES(mathparser ${TARGET_LL})

