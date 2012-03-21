SET(COUPLED_PDE_SRCS
  BasePairCoupling.cc
  CoupledPDEDef.cc
  CouplingMemento.cc
  DirectCoupledPDE.cc
  IterCoupledPDE.cc
  PDECoupling.cc
  PiezoCoupling.cc
  AcouMechCoupling.cc
  )

ADD_LIBRARY(coupledpde STATIC ${COUPLED_PDE_SRCS})

SET(TARGET_LL
  driver
  mathparser
)

TARGET_LINK_LIBRARIES(coupledpde ${TARGET_LL})
