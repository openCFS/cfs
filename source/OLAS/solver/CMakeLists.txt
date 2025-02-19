SET(SOLVER_SRCS
    BaseSolver.cc
    BaseEigenSolver.cc
    generateEigensolver.cc    
    generatesolver.cc
    QuadraticEigenSolver.cc
    ExternalEigenSolver.cc
    )

SET(SOLVER_SRCS
    ${SOLVER_SRCS}
    CGSolver.cc
    DiagSolver.cc
    GMRESSolver.cc
    LDLSolver.cc
    LUSolver.cc
    MINRESSolver.cc
    RichardsonSolver.cc
    ExternalSolver.cc
    )

if(USE_SUPERLU AND USE_ARPACK)
  set(SOLVER_SRCS ${SOLVER_SRCS} PALMSolver.cc)
endif()

ADD_LIBRARY(solver-olas STATIC ${SOLVER_SRCS})

# lapack is always present (netlib, openblas and mkl)
SET(TARGET_LL 
  utils-olas
  lapack-olas )

IF(USE_PARDISO)
  SET(TARGET_LL ${TARGET_LL} pardiso-olas)
ENDIF(USE_PARDISO)

IF(USE_SUITESPARSE)
  SET(TARGET_LL ${TARGET_LL} umfpack-olas)
  INCLUDE_DIRECTORIES(${SUITESPARSE_INCLUDE_DIR})
  SET(TARGET_LL ${TARGET_LL} cholmod-olas)
ENDIF(USE_SUITESPARSE)

IF(USE_SUPERLU)
  INCLUDE_DIRECTORIES(${SUPERLU_INCLUDE_DIR})
  SET(TARGET_LL ${TARGET_LL} superlu-olas)
ENDIF(USE_SUPERLU)

IF(USE_LIS)
  SET(TARGET_LL ${TARGET_LL} lis-olas)
ENDIF(USE_LIS)

if(USE_GINKGO)
  list(APPEND TARGET_LL ginkgo-olas)
endif()

IF(USE_PETSC)
  SET(TARGET_LL ${TARGET_LL} petsc-olas)
ENDIF(USE_PETSC)

IF(USE_ARPACK)
  SET(TARGET_LL ${TARGET_LL} arpack-olas)
ENDIF(USE_ARPACK)

IF(USE_PHIST_EV OR USE_PHIST_CG)
  SET(TARGET_LL ${TARGET_LL} phist-olas  ${PHIST_LIBRARY} ${GHOST_LIBRARY} ${HWLOC_LIBRARY})
  SET_SOURCE_FILES_PROPERTIES(generateEigensolver.cc PROPERTIES COMPILE_FLAGS "-I${MKL_INCLUDE_DIR} ")
ENDIF()


IF(USE_FEAST)
  SET(TARGET_LL ${TARGET_LL} feast-olas)
ENDIF(USE_FEAST)



TARGET_LINK_LIBRARIES(solver-olas ${TARGET_LL})
