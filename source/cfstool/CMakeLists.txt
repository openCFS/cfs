# List of source codes for current target.
SET(CFSTOOL_SRCS cfstool.cc)

# Main CFS++ Tool executable target.
ADD_EXECUTABLE(cfstool ${CFSTOOL_SRCS})


SET(TARGET_LL ${TARGET_LL}
  cmdline
  optimization
  driver
  gridcfs
  pde
  domain
  gridcfs
  pde
  domain
  optimization
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
  SET(TARGET_LL ${TARGET_LL} siminouthdf5)
ENDIF(USE_HDF5)

IF(USE_ANSYSRST)
  SET(TARGET_LL ${TARGET_LL} simoutansysrst)
ENDIF(USE_ANSYSRST)

IF(USE_TCL)
  SET(TARGET_LL ${TARGET_LL} tclmess scripting)
ENDIF(USE_TCL)

IF(USE_PYTHON)
  SET(TARGET_LL ${TARGET_LL} pythonmess scripting)
ENDIF(USE_PYTHON)

# Target link libraries for cfstool target.
SET(TARGET_LL ${TARGET_LL}  
  domain
  elements
  domain
  datainout
  materials
  paramh
  driver
  utils
  matrix
  general
  odesolver
  odedescr
  paramident
  logging
  scripting
  simouttext
  boost_serialization
  boost_filesystem
)
  
SET(TARGET_LL ${TARGET_LL} 
  algsys-olas
  graph-olas
  matvec-olas
  precond-olas
  solver-olas
  utils-olas
)  
  
IF(USE_ARPACK)
  SET(TARGET_LL ${TARGET_LL} arpack-olas)
ENDIF(USE_ARPACK)

IF(USE_ILUPACK)
  SET(TARGET_LL ${TARGET_LL} ilupack-olas)
ENDIF(USE_ILUPACK)

IF(USE_LAPACK)
  SET(TARGET_LL ${TARGET_LL} lapack-olas)
ENDIF(USE_LAPACK)

IF(USE_PARDISO)
  SET(TARGET_LL ${TARGET_LL} pardiso-olas)
ENDIF(USE_PARDISO)
  
TARGET_LINK_LIBRARIES(cfstool ${TARGET_LL})

# Make sure CFS Tool gets linked with system libgcc.so instead of ANSYS libgcc.so.
IF(USE_ANSYSRST)
  SET_TARGET_PROPERTIES(cfstool PROPERTIES LINK_FLAGS -Wl,-rpath,/${LIB_SUFFIX})
ENDIF(USE_ANSYSRST)

# Install the cfstool executable target to ${CMAKE_INSTALL_PREFIX}/bin
INSTALL_TARGETS(/bin cfstool)