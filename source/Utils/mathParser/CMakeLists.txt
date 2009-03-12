# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Attention: Special flag for muParser, as gcc > 4.2 seems to produce faulty
# code with -O3, which could be tracked down to the following compiler flag
#add_definitions(-fno-inline-functions)

SET(MATHPARSER_SRCS mathParser.cc)

ADD_LIBRARY(mathparser STATIC ${MATHPARSER_SRCS})

SET(TARGET_LL
  datainout
  utils
  ${MUPARSER_LIBRARY}
)

TARGET_LINK_LIBRARIES(mathparser ${TARGET_LL})

