SET(OPTIMIZATION_SRCS
  ErsatzMaterial.cc
  Optimization.cc
  Function.cc
  LevelSet.cc
  Condition.cc
  ShapeOpt.cc
  ShapeMapping.cc
  SplineBoxOpt.cc
  TopGrad.cc
  OptimizationMaterial.cc
  Excitation.cc
  StressConstraint.cc
  StateSolution.cc
  SIMP.cc
  PiezoSIMP.cc
  MagSIMP.cc
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
  utils
  algsys-olas
  cfsgeneral
  datainout
  domain
  driver
  logging
  materials
  matvec
  paramh
  pde)

INCLUDE_DIRECTORIES(${XERCES_INCLUDE_DIR})
SET(TARGET_LL ${TARGET_LL} ${XERCES_LIBRARY})

TARGET_LINK_LIBRARIES(optimization ${TARGET_LL})
                   
