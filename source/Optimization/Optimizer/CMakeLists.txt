SET(OPTIMIZER_SRCS
  BaseOptimizer.cc
  ShapeOptimizer.cc
  OptimalityCondition.cc
  EvaluateOnly.cc
  GradientCheck.cc
)

SET(TARGET_LL
  cfsgeneral
  datainout
  design
  domain
  driver
  logging
  matvec
  optimization
  paramh
  pde
  utils
)

IF(USE_IPOPT)
   SET(TARGET_LL ${TARGET_LL} ipopt gfortranbegin gcc_s)

   INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} IPOPT.cc)
ENDIF(USE_IPOPT)

IF(USE_KNITRO)

   #SET(TARGET_LL ${TARGET_LL} ${KNITRO_LIBRARY})
   SET(TARGET_LL ${TARGET_LL} knitro700)
   
   INCLUDE_DIRECTORIES(${KNITRO_INCLUDE_DIR})
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} KNITRO.cc)
ENDIF(USE_KNITRO)

IF(USE_SCPIP)
   SET(TARGET_LL ${TARGET_LL} ${SCPIP_LIBRARY} ${CFS_FORTRAN_LIBS})

   # SCPIPBase is also the open source project C++SCPIP.
   # With this definition the file lives in both worlds! 
   ADD_DEFINITIONS(-DUSE_4_CFS)
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} SCPIP.cc SCPIPBase.cc)
ENDIF(USE_SCPIP)

IF(USE_SNOPT)
   SET(TARGET_LL ${TARGET_LL} ${SNOPT_LIBRARY} ${CFS_FORTRAN_LIBS})
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} SnOptInterface.cc SnOpt.cc)
ENDIF(USE_SNOPT)

ADD_DEFINITIONS(-DUSE_4_CFS)

ADD_LIBRARY(optimizer STATIC ${OPTIMIZER_SRCS})

TARGET_LINK_LIBRARIES(optimizer optimization ${TARGET_LL})
