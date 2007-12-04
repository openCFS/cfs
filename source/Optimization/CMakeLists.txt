SET(OPTIMIZATION_SRCS DesignElement.cc DesignSpace.cc Optimization.cc PiezoSIMP.cc SIMP.cc 
                      TransferFunction.cc Condition.cc OptimalityCondition.cc EvaluateOnly.cc
                      BaseOptimizer.cc)

IF(USE_IPOPT)
   SET(TARGET_LL
    ipopt 
    gfortranbegin 
    gcc_s)

   INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
   
   SET( OPTIMIZATION_SRCS
    ${OPTIMIZATION_SRCS}
    IPOPT.cc)
ENDIF(USE_IPOPT)

IF(USE_SCPIP)
   SET(TARGET_LL
    ${TARGET_LL}
    scpip)

   # SCPIPBase is also the open source project C++SCPIP.
   # With this definition the file lives in both worlds! 
   ADD_DEFINITIONS(-DUSE_4_CFS)
   
   SET( OPTIMIZATION_SRCS
    ${OPTIMIZATION_SRCS}
    SCPIP.cc SCPIPBase.cc)
ENDIF(USE_SCPIP)

ADD_LIBRARY(optimization STATIC ${OPTIMIZATION_SRCS})

SET(TARGET_LL
  ${TARGET_LL}
  domain
  forms
)

TARGET_LINK_LIBRARIES(optimization ${TARGET_LL})

