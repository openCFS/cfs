SET(OPTIMIZATION_SRCS Function.cc Optimization.cc PiezoSIMP.cc ErsatzMaterial.cc SIMP.cc 
                      TransferFunction.cc Condition.cc Objective.cc 
                      ShapeGrad.cc LevelSet.cc ParamMat.cc 
                      TopGrad.cc OptimizationMaterial.cc ShapeOpt.cc )

ADD_DEFINITIONS(-DUSE_4_CFS)

ADD_LIBRARY(optimization STATIC ${OPTIMIZATION_SRCS})

SET(TARGET_LL
  optimizer
  design
  ${TARGET_LL}
  domain
  forms
  matvec
  pde
)

IF(USE_XERCES)
  INCLUDE_DIRECTORIES(${XERCES_INCLUDE_DIR})
  SET(TARGET_LL
    ${TARGET_LL}
    ${XERCES_LIBRARY}
  )
ENDIF(USE_XERCES)


TARGET_LINK_LIBRARIES(optimization ${TARGET_LL})

# --------------------------------------------------------------
# --------- This is for "make scpip" only ----------------------.
# ------------------ Not 4 trunk ! -----------------------------
# SET(SCPIP_SRCS SCPIP_HS071.cc SCPIPBase.cc)
# ADD_EXECUTABLE(scpip ${SCPIP_SRCS})
# INCLUDE_DIRECTORIES(
#  ${INCLUDE_DIRECTORIES} 
#  "/home/fwein/project/cpp_scpip"
#)
#
#TARGET_LINK_LIBRARIES(scpip blas general logging
#                      /space/fwein/packages/scpip/libscpip.a g2c gfortran m)
#
#TARGET_LINK_LIBRARIES(scpip general utils logging scripting utils matvec datainout domain materials utils
#                       lapack blas
#                       /space/fwein/packages/boost_dev/lib/libboost_serialization-gcc.so
#                       /space/fwein/packages/scpip/libscpip.a g2c gfortran m)                      
