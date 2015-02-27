INCLUDE_DIRECTORIES(
  ${HDF5_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR})

IF(USE_CGNS)
  INCLUDE_DIRECTORIES(${CGNS_INCLUDE_DIR})
ENDIF(USE_CGNS)

IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
  LINK_DIRECTORIES(${IPOPT_LINK_DIR})
ENDIF(USE_IPOPT)

# List of source codes for current target.
SET(CFSTOOL_SRCS
  ParamsInit.cc
  HelperFuncs.cc
  WVT.cc
  MatFile.cc
  cfstool.cc
)

# Main CFS++ Tool executable target.
ADD_EXECUTABLE(cfstoolbin ${CFSTOOL_SRCS})

SET(TARGET_LL ${TARGET_LL}
#  optimization
  driver
  mesh 
  pde
  domain
  mesh 
  pde
  domain
#  optimization
  materials
  mathparser
  coupledpde
  febasis
  simouttext
  simoutinfo
)

IF(USE_FFTW)
  LIST(APPEND TARGET_LL ${FFTW_LIBRARY})
ENDIF(USE_FFTW)

IF(USE_GMV)
  LIST(APPEND TARGET_LL siminoutgmv)
ENDIF(USE_GMV)

IF(USE_GIDPOST)
  INCLUDE_DIRECTORIES(${GIDPOST_INCLUDE_DIR})
  LIST(APPEND TARGET_LL simoutgid)
ENDIF(USE_GIDPOST)

IF(USE_MESH)
  LIST(APPEND TARGET_LL siminmesh)
ENDIF(USE_MESH)

IF(USE_UNV)
  LIST(APPEND TARGET_LL simoutunv)
ENDIF(USE_UNV)

IF(USE_HDF5)
  LIST(APPEND TARGET_LL ${HDF5_LT_LIBRARY})
  LIST(APPEND TARGET_LL siminouthdf5)
ENDIF(USE_HDF5)

IF(USE_ANSYSRST)
  LIST(APPEND TARGET_LL simoutansysrst)
ENDIF(USE_ANSYSRST)

# Target link libraries for cfstool target.
LIST(APPEND TARGET_LL
  domain
  mesh
  febasis
  datainout
  materials
  paramh
  driver
  utils
  matvec
  cfsgeneral
#  paramident
  logging
  simouttext
  ${BOOST_SERIALIZATION_LIB}
  ${BOOST_FILESYSTEM_LIB}
  ${BOOST_SYSTEM_LIB}
  )

IF(NOT MINGW)
  LIST(APPEND TARGET_LL
    -lpthread
    )
ENDIF()
  
LIST(APPEND TARGET_LL
  algsys-olas
  graph-olas
  precond-olas
  solver-olas
  utils-olas
)  
  
IF(USE_ARPACK)
  LIST(APPEND TARGET_LL arpack-olas)
ENDIF(USE_ARPACK)

IF(USE_LAPACK)
  LIST(APPEND TARGET_LL lapack-olas)
ENDIF(USE_LAPACK)

IF(USE_PARDISO)
  LIST(APPEND TARGET_LL pardiso-olas)
ENDIF(USE_PARDISO)

LIST(APPEND TARGET_LL
  ${CFS_FORTRAN_DYNRT_LIBS}
  )

TARGET_LINK_LIBRARIES(cfstoolbin ${TARGET_LL})

#IF(MINGW)
#  IF(CFS_BLAS_LAPACK STREQUAL MKL)
#    SET_TARGET_PROPERTIES(cfstoolbin PROPERTIES
#      LINK_FLAGS -Wl,--wrap,__chkstk
#      LINKER_LANGUAGE CXX)
#  ENDIF()
#ENDIF()

SET(CFS_LINK_FLAGS "")

# Make sure CFS Tool gets linked with system libgcc.so instead of ANSYS libgcc.so.
IF(UNIX AND USE_ANSYSRST)
  SET(CFS_LINK_FLAGS "${CFS_LINK_FLAGS} -Wl,-rpath,/${LIB_SUFFIX}")
ENDIF(UNIX AND USE_ANSYSRST)

# Even if we did not compile CFS++ itself with OpenMP enabled, we may have
# compiled external libraries like LIS with OpenMP support.
IF(OPENMP_FOUND)
  SET(CFS_LINK_FLAGS "${CFS_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

SET_TARGET_PROPERTIES(cfstoolbin PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")

# Install the cfstool executable target to ${CMAKE_INSTALL_PREFIX}/bin
#INSTALL_TARGETS(/bin cfstool)
