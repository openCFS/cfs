SET(COUPLED_PDE_SRCS
  IterCoupledPDE.cc
  CoupledPDEDef.cc
  PDECoupling.cc
  CouplingMemento.cc
  DirectCoupledPDE.cc
  BasePairCoupling.cc
  )

ADD_LIBRARY(coupledpde STATIC ${COUPLED_PDE_SRCS})

SET(TARGET_LL
  driver
  mathparser
)

TARGET_LINK_LIBRARIES(coupledpde ${TARGET_LL})
