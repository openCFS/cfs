SET(PRECOND_SRCS
    BasePrecond.cc
    IdPrecondStd.cc
    generateprecond.cc
    )

SET(PRECOND_SRCS
    ${PRECOND_SRCS}
    IC0Precond.cc
    ILU0Precond.cc
    ILUK_Precond.cc
    ILUTP_Precond.cc
    JacPrecond.cc
    MGPrecond.cc
    SSORPrecond.cc
    SBMDiagPrecond.cc
    SBMJacobiPrecond.cc
    ILDLPrecond/ildl0factoriser.cc
    ILDLPrecond/ildlcnfactoriser.cc
    ILDLPrecond/ildlkfactoriser.cc
    ILDLPrecond/ildlprecond.cc
    ILDLPrecond/ildltpfactoriser.cc
    LDLPrecond.cc
    )

ADD_LIBRARY(precond-olas STATIC ${PRECOND_SRCS})

SET(TARGET_LL
    matvec
    utils-olas
    )
  
IF(MKL_INCLUDE_DIR)
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
ENDIF(MKL_INCLUDE_DIR)


TARGET_LINK_LIBRARIES(precond-olas ${TARGET_LL})
