SET(MATVEC_SRCS
  BaseMatrix.cc
  generatematvec.cc
  SBM_Vector.cc
  SBM_Matrix.cc
  StdMatrix.cc
  matrixSerialization.cc
  vectorSerialization.cc
  CRS_Matrix.cc
  Diag_Matrix.cc
  SCRS_Matrix.cc
  Vector.cc
  VBR_Matrix.cc
  Matrix.cc
)

IF(MKL_INCLUDE_DIR)
  INCLUDE_DIRECTORIES(SYSTEM ${MKL_INCLUDE_DIR})
ENDIF(MKL_INCLUDE_DIR)

ADD_LIBRARY(matvec STATIC ${MATVEC_SRCS})

SET(TARGET_LL
  ${TARGET_LL}
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS})

SET(TARGET_LL
  ${TARGET_LL}
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS})

SET(TARGET_LL
  ${TARGET_LL}
  ${BOOST_SERIALIZATION_LIB}
  graph-olas)

TARGET_LINK_LIBRARIES(matvec ${TARGET_LL})
