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
  CoordFormat.cc
)

IF(MKL_INCLUDE_DIR)
  INCLUDE_DIRECTORIES(SYSTEM ${MKL_INCLUDE_DIR})
ENDIF(MKL_INCLUDE_DIR)

if(USE_EMBEDDED_PYTHON)
  # see FindPrograms.cmake on how the PYTHON_ stuff is set
  list(APPEND MATVEC_SRCS VectorPython.cc)
  list(APPEND TARGET_LL ${Python_LIBRARIES})
endif()

ADD_LIBRARY(matvec STATIC ${MATVEC_SRCS})

SET(TARGET_LL
  ${TARGET_LL}
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS}
  ${BOOST_LIBRARY}
  ${ZLIB_LIBRARY} 
  graph-olas)

TARGET_LINK_LIBRARIES(matvec ${TARGET_LL})
