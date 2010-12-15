SET(MATVEC_SRCS
    basematrix.cc
    generatematvec.cc
    sbmvector.cc
    sbmmatrix.cc
    SingleVector.cc
    stdmatrix.cc
    matrixSerialization.cc
    vectorSerialization.cc
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

IF(USE_LAPACK)
  SET(TARGET_LL
  ${TARGET_LL}
  lapack-olas
  )
ENDIF(USE_LAPACK)

SET(TARGET_LL
    ${TARGET_LL}
    ${BOOST_SERIALIZATION_LIB}
)

TARGET_LINK_LIBRARIES(matvec ${TARGET_LL})
