IF(USE_MPCCI)
  INCLUDE_DIRECTORIES(${MPCCI_INCLUDE_DIR})
ENDIF(USE_MPCCI)

IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
  LINK_DIRECTORIES(${IPOPT_LINK_DIR})
ENDIF(USE_IPOPT)

# List of source codes for current target.
SET(CFS_SRCS cfs.cc)

# Main CFS++ executable target.
ADD_EXECUTABLE(cfs ${CFS_SRCS})

# Define a list of libraries CFS++ depends upon. At the moment there
# are some nasty circular dependencies which should be resolved in
# the future (see e.g. pde <-> domain, domain <-> gridcfs !!!)
  
#IF(USE_MPCCI)
#  SET(TARGET_LL ${TARGET_LL} mpccicpl)
#ENDIF(USE_MPCCI)

SET(TARGET_LL
  general
  utils
  cmdline
  definefiles
  paramh
  domain
  )

# Target link libraries for CFS++ executable target.
TARGET_LINK_LIBRARIES(cfs 
  ${TARGET_LL}
  )


# Make sure CFS++ gets linked with system libgcc.so instead of ANSYS libgcc.so.
IF(USE_ANSYSRST)
  SET_TARGET_PROPERTIES(cfs PROPERTIES LINK_FLAGS -Wl,-rpath,/${LIB_SUFFIX})
ENDIF(USE_ANSYSRST)

# Install the CFS++ executable target to ${CMAKE_INSTALL_PREFIX}/bin
INSTALL_TARGETS(/bin cfs)
