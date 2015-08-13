INCLUDE_DIRECTORIES(
  ${HDF5_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR})

# List of source codes for current target.
SET(CFSDAT_SRCS
  cfsdat.cc
  Filters/Input/InputFilterGeneric.cc
)

# Main CFS++ Dat executable target.
ADD_EXECUTABLE(cfsdatbin ${CFSDAT_SRCS})

SET(TARGET_LL 
 datainout
)


TARGET_LINK_LIBRARIES(cfsdatbin ${TARGET_LL})

SET(CFS_LINK_FLAGS "")

SET_TARGET_PROPERTIES(cfsdatbin PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")

# Install the cfstool executable target to ${CMAKE_INSTALL_PREFIX}/bin
#INSTALL_TARGETS(/bin cfsdat)