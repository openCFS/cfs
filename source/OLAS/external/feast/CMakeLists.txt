SET(FEAST_SRCS FeastEigenSolver.cc)

ADD_LIBRARY(feast-olas STATIC ${FEAST_SRCS})
# as long as not the academic FEAST implementation is used we need MKL which contains FEAST
SET(TARGET_LL ${FEAST_LIBRARY} ${BLAS_LIBRARY} ${CFS_FORTRAN_LIBS})

TARGET_LINK_LIBRARIES(feast-olas ${TARGET_LL})

ADD_DEPENDENCIES(feast-olas feast boost)

INCLUDE_DIRECTORIES(${FEAST_INCLUDE_DIR})
