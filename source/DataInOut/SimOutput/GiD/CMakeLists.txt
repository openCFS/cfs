# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(SIMOUTGID_SRCS simOutGiD.cc)

ADD_LIBRARY(simoutgid STATIC ${SIMOUTGID_SRCS})

TARGET_LINK_LIBRARIES(simoutgid gidpost)
