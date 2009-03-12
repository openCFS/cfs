# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(SIMOUTUNV_SRCS
  simOutUnv.cc
  simInputUnv.cc
  unv_reader.cc
  unv.cc
  unv_dat.cc
  unv_if.cc
)

ADD_LIBRARY(simoutunv STATIC ${SIMOUTUNV_SRCS})
