SET(OPTIMIZER_SRCS
  BaseOptimizer.cc
  ShapeOptimizer.cc
  OptimalityCondition.cc
  EvaluateOnly.cc
  GradientCheck.cc
  SGP.cc
)

SET(TARGET_LL
  optimization
  cfsgeneral
  datainout
  design
  domain
  driver
  logging
  matvec
  paramh
  pde
  utils
)

# we need it almost always ?!
SET(DEPENDS
  "${LAPACK_LIBRARY}"
  "${BLAS_LIBRARY}"
  "${CFS_FORTRAN_LIBS}"
)

IF(USE_KNITRO)

   #SET(TARGET_LL ${TARGET_LL} ${KNITRO_LIBRARY})
   SET(TARGET_LL ${TARGET_LL} knitro700)
   
   INCLUDE_DIRECTORIES(${KNITRO_INCLUDE_DIR})
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} KNITRO.cc)
ENDIF(USE_KNITRO)

IF(USE_SCPIP)
   SET(TARGET_LL ${TARGET_LL} ${SCPIP_LIBRARY} ${DEPENDS})
   # message("${CFS_FORTRAN_LIBS}")
   # SCPIPBase is also the open source project C++SCPIP.
   # With this definition the file lives in both worlds! 
   ADD_DEFINITIONS(-DUSE_4_CFS)
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} SCPIP.cc SCPIPBase.cc)
ENDIF(USE_SCPIP)

ADD_DEFINITIONS(-DUSE_4_CFS)

IF(USE_SNOPT)
   SET(TARGET_LL ${TARGET_LL}  ${SNOPT_LIBRARY} ${DEPENDS})
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} SnOpt.cc)
   SET(OPTIMIZER_FSRCS SnOptInterface.f)
   ADD_LIBRARY(snopt_ftnif STATIC ${OPTIMIZER_FSRCS})
ENDIF(USE_SNOPT)

IF(USE_IPOPT)
   SET(TARGET_LL ${TARGET_LL}  ${IPOPT_LIBRARY} ${METIS_LIBRARY} ${DEPENDS})
   INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} IPOPT.cc FeasPP.cc FeasSubProblem.cc)
ENDIF(USE_IPOPT)


ADD_LIBRARY(optimizer STATIC ${OPTIMIZER_SRCS})

IF(USE_SNOPT)
  TARGET_LINK_LIBRARIES(optimizer snopt_ftnif ${TARGET_LL})
ELSE()
  TARGET_LINK_LIBRARIES(optimizer ${TARGET_LL})
ENDIF()
