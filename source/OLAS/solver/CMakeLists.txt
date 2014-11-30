# Make sure the compiler can find include files
INCLUDE_DIRECTORIES(${CFSDEPS_INCLUDE_DIR})

SET(SOLVER_SRCS
  basesolver.cc
  generatesolver.cc
  baseEigensolver.cc
  generateEigensolver.cc
    )

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(SOLVER_SRCS
    ${SOLVER_SRCS}
    cgsolver.cc
    diagsolver.cc
    gmres.cc
    ldlsolver.cc
    lusolver.cc
    minres.cc
    richardson.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(solver-olas STATIC ${SOLVER_SRCS})

SET(TARGET_LL
  cfsgeneral
  lapack-olas
  matvec
  paramh
  precond-olas
  utils
  utils-olas
)

IF(USE_PARDISO)
  SET(TARGET_LL ${TARGET_LL} pardiso-olas)
ENDIF(USE_PARDISO)

IF(USE_PARDISO64)
  SET(TARGET_LL ${TARGET_LL} pardiso64-olas)
ENDIF(USE_PARDISO64)

IF(USE_ILUPACK)
  SET(TARGET_LL ${TARGET_LL} ilupack-olas)
ENDIF(USE_ILUPACK)

IF(USE_CHOLMOD)
  INCLUDE_DIRECTORIES(${CHOLMOD_INCLUDE_DIR})
  SET(TARGET_LL ${TARGET_LL} cholmod-olas)
ENDIF(USE_CHOLMOD)

IF(USE_IPOPT)
  SET(TARGET_LL ${TARGET_LL} contact-ipopt)
ENDIF(USE_IPOPT)

IF(USE_ARPACK)
  SET(TARGET_LL ${TARGET_LL} arpack-olas)
ENDIF(USE_ARPACK)


TARGET_LINK_LIBRARIES(solver-olas ${TARGET_LL})
