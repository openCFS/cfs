SET(PRECOND_SRCS
    baseprecond.cc
    idprecond.cc
    generateprecond.cc
    #mgprecond.cc
    #sbmdiagprecond.cc
    )

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(PRECOND_SRCS
    ${PRECOND_SRCS}
    ic0precond.cc
    ilu0precond.cc
    ilukprecond.cc
    ilutpprecond.cc
    jacprecond.cc
    ssorprecond.cc
    ILDLPrecond/ildl0factoriser.cc
    ILDLPrecond/ildlcnfactoriser.cc
    ILDLPrecond/ildlkfactoriser.cc
    ILDLPrecond/ildlprecond.cc
    ILDLPrecond/ildltpfactoriser.cc    
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(precond-olas STATIC ${PRECOND_SRCS})

SET(TARGET_LL
    matvec
    utils-olas
    )
    
TARGET_LINK_LIBRARIES(precond-olas ${TARGET_LL})
