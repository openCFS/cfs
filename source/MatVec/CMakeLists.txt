SET(MATVEC_SRCS
    BaseMatrix.cc
    generatematvec.cc
    SBM_Vector.cc
    SBM_Matrix.cc
    SingleVector.cc
    StdMatrix.cc
    matrixSerialization.cc
    vectorSerialization.cc
    )

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(MATVEC_SRCS
    ${MATVEC_SRCS}
    CRS_Matrix.cc
    Diag_Matrix.cc
    SCRS_Matrix.cc
    TypeDefs.cc
    Vector.cc
    VBR_Matrix.cc
    Matrix.cc
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
    graph-olas
)

TARGET_LINK_LIBRARIES(matvec ${TARGET_LL})
