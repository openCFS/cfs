SET(COUPLED_PDE_SRCS
  BasePairCoupling.cc
  DirectCoupledPDE.cc
  IterCoupledPDE.cc
  PiezoCoupling.cc
  MagnetoStrictCoupling.cc
  AcouMechCoupling.cc
  AcouMechCoupling.cc
  FluidMechCoupling.cc
  WaterWaveAcousticsCoupling.cc
  WaterWaveMechCoupling.cc
  LinFlowHeatCoupling.cc
  LinFlowAcouCoupling.cc
  LinFlowMechCoupling.cc
  )

ADD_LIBRARY(coupledpde STATIC ${COUPLED_PDE_SRCS})

SET(TARGET_LL
  driver
  mathparser
)

TARGET_LINK_LIBRARIES(coupledpde ${TARGET_LL})
