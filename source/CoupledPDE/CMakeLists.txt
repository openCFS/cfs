SET(COUPLED_PDE_SRCS
  PiezoCoupling.cc
  BasePairCoupling.cc
  DirectCoupledPDE.cc
  coupledpdedef.cc
  pdecoupling.cc
  itercoupledpde.cc
  magStrictCoupling.cc
  AcouMechCoupling.cc
  ThermoMechCoupling.cc
  ThermoElectricCoupling.cc
  couplingmemento.cc
  )

ADD_LIBRARY(coupledpde STATIC ${COUPLED_PDE_SRCS})

TARGET_LINK_LIBRARIES(coupledpde
  algsys-olas
  cfsgeneral
  datainout
  domain
  driver
  elements
  forms
  logging
  materials
  matvec
  paramh
  pde
  utils  
)
