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

# Main cfstool executable target.
ADD_EXECUTABLE(cfstool ${CFSTOOL_SRCS})

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

IF(USE_GIDPOST)
  INCLUDE_DIRECTORIES(${GIDPOST_INCLUDE_DIR})
  LIST(APPEND TARGET_LL simoutgid)
ENDIF(USE_GIDPOST)

LIST(APPEND TARGET_LL simoutunv)

# Target link libraries for cfstool target.
# add optimization to avoid a strange linke problem introduced in @13977 
LIST(APPEND TARGET_LL
  ${HDF5_LIBRARY}
  siminouthdf5
  siminmesh
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
  ${BOOST_LIBRARY}
  ${ZLIB_LIBRARY} )

IF(NOT WIN32)
  LIST(APPEND TARGET_LL -lpthread ${CFS_FORTRAN_DYNRT_LIBS} )
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

if(USE_BLAS_LAPACK STREQUAL "NETLIB")
  assert_set(LAPACK_LIBRARY)
  list(APPEND TARGET_LL ${LAPACK_LIBRARY})
  list(APPEND TARGET_LL ${BLAS_LIBRARY})
elseif(USE_BLAS_LAPACK STREQUAL "MKL")
  list(APPEND TARGET_LL ${MKL_BLAS_LIB})
endif()

LIST(APPEND TARGET_LL lapack-olas)

IF(USE_PARDISO)
  LIST(APPEND TARGET_LL pardiso-olas)
ENDIF(USE_PARDISO)

LIST(APPEND TARGET_LL ${CFS_FORTRAN_DYNRT_LIBS} )

TARGET_LINK_LIBRARIES(cfstool ${TARGET_LL})

IF(CMAKE_GENERATOR STREQUAL "Visual Studio 16 2019")
  SET(INTEL_LIB_DIR "${INTEL_COMPILER_DIR}/../../compiler/lib/intel64_win")
  TARGET_LINK_DIRECTORIES(cfstool PUBLIC ${INTEL_LIB_DIR})
ENDIF()

SET(CFS_LINK_FLAGS "")

# Even if we did not compile openCFS itself with OpenMP enabled, we may have
# compiled external libraries like LIS with OpenMP support.
IF(OPENMP_FOUND)
  SET(CFS_LINK_FLAGS "${CFS_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

SET_TARGET_PROPERTIES(cfstool PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")

# IntelLLVM und RHEL8 chooses a too old lib for std::filesystem. CFS_LINKER_FLAGS does not work
if(UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  target_link_libraries(cfstool stdc++fs)
endif()



