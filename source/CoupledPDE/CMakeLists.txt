SET(COUPLED_PDE_SRCS itercoupledpde.cc coupledpdedef.cc pdecoupling.cc
  couplingmemento.cc DirectCoupledPDE.cc BasePairCoupling.cc
  PiezoCoupling.cc AcouMechCoupling.cc ThermoMechCoupling.cc
  ThermoElectricCoupling.cc)

ADD_LIBRARY(coupledpde STATIC ${COUPLED_PDE_SRCS})

SET(TARGET_LL
  driver
  forms
)

TARGET_LINK_LIBRARIES(coupledpde ${TARGET_LL})
