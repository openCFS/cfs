# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(SIMOUTTEXT_SRCS textSimOutput.cc)

ADD_LIBRARY(simouttext STATIC ${SIMOUTTEXT_SRCS})

SET(TARGET_LL
  datainout
)

TARGET_LINK_LIBRARIES(simouttext ${TARGET_LL})
