SET(FEAST_SRCS FeastEigenSolver.cc)

ADD_LIBRARY(feast-olas STATIC ${FEAST_SRCS})

# as long as not the academic FEAST implementation is used we need MKL which contains FEAST
SET(TARGET_LL ${FEAST_LIBRARY} ${BLAS_LIBRARY} ${CFS_FORTRAN_LIBS})

# fix for MKL >= 2026, source from Eric Polizzi, should not be needed any more once FEAST is updated > 4.0.
if(USE_BLAS_LAPACK STREQUAL "MKL")
  include_directories(${MKL_INCLUDE_DIR})
  add_library(feast-c STATIC "${CMAKE_SOURCE_DIR}/cfsdeps/feast/wrapper_mkl_sparse.c")
  list(APPEND TARGET_LL "feast-c")
endif()

TARGET_LINK_LIBRARIES(feast-olas ${TARGET_LL})

ADD_DEPENDENCIES(feast-olas feast boost)

INCLUDE_DIRECTORIES(${FEAST_INCLUDE_DIR})
