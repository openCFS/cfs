# Make sure the compiler can find include files from our CFSMESHIO library.
INCLUDE_DIRECTORIES(${MPCCI_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
# LINK_DIRECTORIES()

SET(MPCCICPL_SRCS MpCCIexch.cc)

ADD_LIBRARY(mpccicpl STATIC ${MPCCICPL_SRCS})

TARGET_LINK_LIBRARIES(mpccicpl ${MPCCI_LIBRARIES})