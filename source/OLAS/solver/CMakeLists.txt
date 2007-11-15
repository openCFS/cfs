# Make sure the compiler can find include files
INCLUDE_DIRECTORIES(${CFSDEPS_INCLUDE_DIR})

SET(SOLVER_SRCS generatesolver.cc generateEigensolver.cc)

ADD_DEFINITIONS()

ADD_LIBRARY(solver-olas STATIC ${SOLVER_SRCS})


SET(TARGET_LL
  utils-olas
)

IF(USE_PARDISO)
  SET(TARGET_LL ${TARGET_LL} pardiso-olas)
ENDIF(USE_PARDISO)

IF(USE_ILUPACK)
  SET(TARGET_LL ${TARGET_LL} ilupack-olas)
ENDIF(USE_ILUPACK)

IF(USE_ARPACK)
  SET(TARGET_LL ${TARGET_LL} arpack-olas)
ENDIF(USE_ARPACK)


TARGET_LINK_LIBRARIES(solver-olas ${TARGET_LL})
