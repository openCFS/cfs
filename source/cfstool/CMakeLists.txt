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
  driver
  mesh 
  pde
  domain
  mesh 
  pde
  domain
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

# Target link libraries for cfstool target.
# add optimization to avoid a strange linke problem introduced in @13977 
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
  optimization
  logging
  simouttext
  ${BOOST_SERIALIZATION_LIB}
  ${BOOST_FILESYSTEM_LIB}
  )

IF(NOT WIN32)
  IF(NOT MINGW)
    LIST(APPEND TARGET_LL -lpthread ${CFS_FORTRAN_DYNRT_LIBS} )
  ENDIF()
ENDIF(NOT WIN32)
  
LIST(APPEND TARGET_LL
  algsys-olas
  graph-olas
  precond-olas
  solver-olas
  utils-olas )  
  
IF(USE_ARPACK)
  LIST(APPEND TARGET_LL arpack-olas)
ENDIF(USE_ARPACK)

IF(USE_PHIST_EV OR USE_PHIST_CG)
  LIST(APPEND TARGET_LL ${CUDA_LIBS} ${HWLOC_LIBRARY} )
ENDIF(USE_PHIST_EV OR USE_PHIST_CG )

IF(USE_ILUPACK OR USE_ILUPACK_PARALLEL)
  LIST(APPEND TARGET_LL ${ILUPACK_LIBRARY} )
ENDIF(USE_ILUPACK OR USE_ILUPACK_PARALLEL )


LIST(APPEND TARGET_LL lapack-olas)

IF(USE_PARDISO)
  LIST(APPEND TARGET_LL pardiso-olas)
ENDIF(USE_PARDISO)

#IF(MKL_BLAS_LIB)
#  LIST(APPEND TARGET_LL ${MKL_BLAS_LIB})
#ENDIF(MKL_BLAS_LIB)

##IF(DEPS_SEQUENTIAL)
#	#LIST(APPEND TARGET_LL ${DEPS_SEQUENTIAL})
#  #ENDIF(DEPS_SEQUENTIAL)

LIST(APPEND TARGET_LL ${CFS_FORTRAN_DYNRT_LIBS} )

TARGET_LINK_LIBRARIES(cfstoolbin ${TARGET_LL})

#IF(MINGW)
#  IF(CFS_BLAS_LAPACK STREQUAL MKL)
#    SET_TARGET_PROPERTIES(cfstoolbin PROPERTIES
#      LINK_FLAGS -Wl,--wrap,__chkstk
#      LINKER_LANGUAGE CXX)
#  ENDIF()
#ENDIF()

SET(CFS_LINK_FLAGS "")

# Even if we did not compile CFS++ itself with OpenMP enabled, we may have
# compiled external libraries like LIS with OpenMP support.
IF(OPENMP_FOUND)
  SET(CFS_LINK_FLAGS "${CFS_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

SET_TARGET_PROPERTIES(cfstoolbin PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")
