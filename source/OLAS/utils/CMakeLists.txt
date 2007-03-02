# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(UTILS_SRCS environment-olas.cc output.cc tools-olas.cc profiler-olas.cc
  Exception.cc Enum.cc math/givensrotation.cc math/iterativerefinement.cc)

ADD_LIBRARY(utils-olas STATIC ${UTILS_SRCS})


