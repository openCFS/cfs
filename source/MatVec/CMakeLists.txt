SET(MATVEC_SRCS
  vectorSerialization.cc
  generatematvec.cc
  stdmatrix.cc
  sbmmatrix.cc
  SingleVector.cc
  sbmvector.cc
  matrixSerialization.cc
  basematrix.cc
    )

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(MATVEC_SRCS
    ${MATVEC_SRCS}
    crs_matrix.cc
    diag_matrix.cc
    scrs_matrix.cc
    typedefs.cc
    vector.cc
    matrix.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(matvec STATIC ${MATVEC_SRCS})

SET(TARGET_LL
  ${BOOST_SERIALIZATION_LIB}
  cfsgeneral
  graph-olas
  lapack-olas
  utils  
)
IF(USE_LAPACK)
  SET(TARGET_LL ${TARGET_LL} lapack-olas )
ENDIF(USE_LAPACK)

TARGET_LINK_LIBRARIES(matvec ${TARGET_LL})
