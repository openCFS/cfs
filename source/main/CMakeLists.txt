
IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
  LINK_DIRECTORIES(${IPOPT_LINK_DIR})
ENDIF(USE_IPOPT)

# List of source codes for current target.
SET(CFS_SRCS CFS.cc)

# Main CFS++ executable target.
ADD_EXECUTABLE(cfsbin ${CFS_SRCS})

# Define a list of libraries CFS++ depends upon. At the moment there
# are some nasty circular dependencies which should be resolved in
# the future (see e.g. pde <-> domain, domain <-> gridcfs !!!)
  
SET(TARGET_LL
  cfsgeneral
  utils
  datainout
  paramh
  utils
  domain
  elemmapping
  ${CFS_FORTRAN_DYNRT_LIBS} )

IF(USE_XERCES)
  LIST(APPEND TARGET_LL ${XERCES_LIBRARY})
ENDIF(USE_XERCES)

IF(USE_LIBXML2)
  LIST(APPEND TARGET_LL ${LIBXML2_LIBRARY})
ENDIF(USE_LIBXML2)  

IF(USE_LIBFBI)
  LIST(APPEND TARGET_LL ${BOOST_THREAD_LIB} )
ENDIF()

IF(USE_PHIST_EV OR USE_PHIST_CG)
  LIST(APPEND TARGET_LL ${CUDA_LIBS} ${HWLOC_LIBRARY})
ENDIF()

IF(MKL_BLAS_LIB)
  LIST(APPEND TARGET_LL ${MKL_BLAS_LIB})
ENDIF(MKL_BLAS_LIB)

IF(DEPS_SEQUENTIAL)
  LIST(APPEND TARGET_LL ${DEPS_SEQUENTIAL})
ENDIF(DEPS_SEQUENTIAL)

# lapack is in NETLIB, OPENBLAS and MKL
# LIST(APPEND TARGET_LL -lpthread -lquadmath )

SET(CFS_LINK_FLAGS "")

# Even if we did not compile CFS++ itself with OpenMP enabled, we may have
# compiled external libraries like LIS with OpenMP support.
IF(OPENMP_FOUND)
  SET(CFS_LINK_FLAGS "${CFS_LINK_FLAGS} ${OpenMP_CXX_FLAGS}  ")
ENDIF()

IF(CMAKE_GENERATOR STREQUAL "Visual Studio 16 2019")
  SET(INTEL_LIB_DIR "${INTEL_COMPILER_DIR}/../../compiler/lib/intel64_win")
  TARGET_LINK_DIRECTORIES(cfsbin PUBLIC ${INTEL_LIB_DIR})
ENDIF()

# Target link libraries for CFS++ executable target.
TARGET_LINK_LIBRARIES(cfsbin ${TARGET_LL})

SET_TARGET_PROPERTIES(cfsbin PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")
