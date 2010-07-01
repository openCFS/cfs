SET(OPTIMIZER_SRCS ShapeOptimizer.cc BaseOptimizer.cc GradientCheck.cc 
                   OptimalityCondition.cc EvaluateOnly.cc)
IF(USE_IPOPT)
   SET(TARGET_LL ipopt gfortranbegin gcc_s)

   INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} IPOPT.cc)
ENDIF(USE_IPOPT)

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
