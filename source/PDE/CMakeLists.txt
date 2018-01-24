IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

SET(PDE_SRCS
  BasePDE.cc
  ElecPDE.cc
  PerturbedFlowPDE.cc
  HeatPDE.cc
  MagEdgePDE.cc
  MagneticPDE.cc
  MechPDE.cc
  SinglePDE.cc
  StdPDE.cc
  AcousticPDE.cc
  AcousticSplitPDE.cc
  AcousticMixedPDE.cc
  ElecCurrentPDE.cc
  WaterWavePDE.cc
  FlowPDE.cc
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
  
)

IF(USE_HDF5)
  SET(TARGET_LL
    ${TARGET_LL}
    siminouthdf5
    )
ENDIF(USE_HDF5)

TARGET_LINK_LIBRARIES(pde ${TARGET_LL})
IF(USE_ADOLC)
  INCLUDE_DIRECTORIES(${ADOLC_INCLUDE_DIR})
  #LINK_DIRECTORIES(${ADOLC_LIBRARY})
ENDIF(USE_ADOLC)


