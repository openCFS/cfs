SET(OPTIMIZATION_SRCS
  ErsatzMaterial.cc
  Optimization.cc
  Function.cc
  LevelSet.cc
  Condition.cc
  ShapeOpt.cc
  ShapeMapping.cc
  TopGrad.cc
  OptimizationMaterial.cc
  Excitation.cc
  StressConstraint.cc
  StateSolution.cc
  SIMP.cc
  PiezoSIMP.cc
  LBMSIMP.cc
  PiezoParamMat.cc
  ShapeGrad.cc
  TransferFunction.cc
  Transform.cc
  Objective.cc
  Context.cc
  ParamMat.cc
 )

ADD_DEFINITIONS(-DUSE_4_CFS)

ADD_LIBRARY(optimization STATIC ${OPTIMIZATION_SRCS})

SET(TARGET_LL
  algsys-olas
  cfsgeneral
  datainout
  domain
  driver
  logging
  materials
  matvec
  paramh
  pde
  utils
)

INCLUDE_DIRECTORIES(${XERCES_INCLUDE_DIR})
SET(TARGET_LL ${TARGET_LL} ${XERCES_LIBRARY})

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
