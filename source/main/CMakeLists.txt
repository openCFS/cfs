IF(USE_MPCCI)
  INCLUDE_DIRECTORIES(${MPCCI_INCLUDE_DIR})
ENDIF(USE_MPCCI)

IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
  LINK_DIRECTORIES(${IPOPT_LINK_DIR})
ENDIF(USE_IPOPT)

# List of source codes for current target.
SET(CFS_SRCS cfs.cc)

IF(DEBUG)
  SET(CFS_EXE cfs)
ELSE(DEBUG)
  SET(CFS_EXE cfs_rel)
ENDIF(DEBUG)

# Main CFS++ executable target.
ADD_EXECUTABLE(${CFS_EXE} ${CFS_SRCS})

# Define a list of libraries CFS++ depends upon. At the moment there
# are some nasty circular dependencies which should be resolved in
# the future (see e.g. pde <-> domain, domain <-> gridcfs !!!)
  
IF(USE_MPCCI)
  SET(TARGET_LL ${TARGET_LL} mpccicpl)
ENDIF(USE_MPCCI)


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
  SET(TARGET_LL ${TARGET_LL} siminoutxmdf)
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

SET(TARGET_LL ${TARGET_LL} 
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

# Target link libraries for CFS++ executable target.
TARGET_LINK_LIBRARIES(${CFS_EXE} ${TARGET_LL}
  boost_serialization
  boost_filesystem)

# Make sure CFS++ gets linked with system libgcc.so instead of ANSYS libgcc.so.
IF(USE_ANSYSRST)
  SET_TARGET_PROPERTIES(cfs PROPERTIES LINK_FLAGS -Wl,-rpath,/${LIB_SUFFIX})
ENDIF(USE_ANSYSRST)

# Install the CFS++ executable target to ${CMAKE_INSTALL_PREFIX}/bin
INSTALL_TARGETS(/bin ${CFS_EXE})