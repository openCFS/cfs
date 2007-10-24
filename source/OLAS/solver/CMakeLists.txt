# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

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
