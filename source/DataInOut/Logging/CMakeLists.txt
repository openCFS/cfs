# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(LOGGING_SRCS functions.cc log.cc log_manager.cc cfslog.cc)

ADD_LIBRARY(logging STATIC ${LOGGING_SRCS})
