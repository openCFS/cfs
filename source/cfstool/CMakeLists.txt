# List of source codes for current target.
SET(CFSTOOL_SRCS cfstool.cc)

# Main CFS++ Tool executable target.
ADD_EXECUTABLE(cfstool ${CFSTOOL_SRCS})

SET(TARGET_LL ${TARGET_LL}
  gridcfs
  pde
  domain
  gridcfs
  pde
  domain
  materials
  mathparser
  coupledpde
  forms
  matrix
  elements
  definefiles
  simouttext
)

IF(USE_GMV)
  IF(USE_GMV_INPUT)
    SET(TARGET_LL ${TARGET_LL} simingmv)
  ENDIF(USE_GMV_INPUT)

  IF(USE_GMV_OUTPUT)
    SET(TARGET_LL ${TARGET_LL} simoutgmv)
  ENDIF(USE_GMV_OUTPUT)
ENDIF(USE_GMV)

IF(USE_GIDPOST)
  SET(TARGET_LL ${TARGET_LL} simoutgid)
ENDIF(USE_GIDPOST)

IF(USE_MESH)
  SET(TARGET_LL ${TARGET_LL} siminmesh)
ENDIF(USE_MESH)

IF(USE_UNV)
  SET(TARGET_LL ${TARGET_LL} simoutunv)
ENDIF(USE_UNV)

IF(USE_HDF5)
  SET(TARGET_LL ${TARGET_LL} siminoutxmdf)
ENDIF(USE_HDF5)

IF(USE_ANSYSRST)
  SET(TARGET_LL ${TARGET_LL} simoutansysrst)
ENDIF(USE_ANSYSRST)

# Target link libraries for CFS++ executable target.
TARGET_LINK_LIBRARIES(cfstool ${TARGET_LL}
  boost_serialization
  boost_filesystem)

# Make sure CFS Tool gets linked with system libgcc.so instead of ANSYS libgcc.so.
IF(USE_ANSYSRST)
  SET_TARGET_PROPERTIES(cfstool PROPERTIES LINK_FLAGS -Wl,-rpath,/${LIB_SUFFIX})
ENDIF(USE_ANSYSRST)

# Install the cfstool executable target to ${CMAKE_INSTALL_PREFIX}/bin
INSTALL_TARGETS(/bin cfstool)