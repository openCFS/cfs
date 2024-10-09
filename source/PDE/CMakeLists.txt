IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

SET(PDE_SRCS
  BasePDE.cc
  ElecPDE.cc
  MagneticScalarPotentialPDE.cc
  MagneticScalarPotentialAdjPDE.cc
  PerturbedFlowPDE.cc
  HeatPDE.cc
  MagBasePDE.cc
  MagEdgePDE.cc
  MagEdgeAdjPDE.cc
  MagEdgeSpecialAVPDE.cc
  MagEdgeMixedAVPDE.cc
  MagEdgeMixedSFGPDE.cc
  MagEdgeRegulSFGPDE.cc
  MagneticPDE.cc
  DarwinPDE.cc
  MechPDE.cc
  SmoothPDE.cc
  SinglePDE.cc
  StdPDE.cc
  AcousticPDE.cc
  AcousticSplitPDE.cc
  AcousticMixedPDE.cc
  ElecCurrentPDE.cc
  ElecQuasiStaticPDE.cc
  WaterWavePDE.cc
  FlowPDE.cc
  LinFlowPDE.cc
  LinFlowPDE.cc
  TestPDE.cc
  LatticeBoltzmannPDE.cc
)

ADD_LIBRARY(pde STATIC ${PDE_SRCS})
INCLUDE_DIRECTORIES(${LBMSOLVER_DIR})

SET(TARGET_LL
  lbm_solver
  datainout
  domain
  algsys-olas
  driver
  mathparser
  paramh
  utils  
  matvec
  timeschemes
  mesh
  results
  febasis
  coeffunction
  linforms
  bilinforms
  siminouthdf5
)

TARGET_LINK_LIBRARIES(pde ${TARGET_LL})
