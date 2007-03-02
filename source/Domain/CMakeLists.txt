# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(DOMAIN_SRCS domain.cc grid.cc entityList.cc resultInfo.cc
  ansatzFct.cc bcs.cc edgeFace.cc)

ADD_LIBRARY(domain STATIC ${DOMAIN_SRCS})
