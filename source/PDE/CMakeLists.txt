IF(USE_MPCCI)
  INCLUDE_DIRECTORIES(${MPCCI_INCLUDE_DIR})
ENDIF(USE_MPCCI)

IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

SET(PDE_SRCS
  mechPDE.cc
  acousticPDE.cc
  magneticPDE.cc
  elecPDE.cc
  fluidMechPDE.cc
  eqnMap.cc
  magEdgePDE.cc
  acousticMixedPDE.cc
  disContEqnMap.cc
  heatCondPDE.cc
  StdPDE.cc
  pdememento.cc
  magneticScalarPDE.cc
  newmarkFracDampMech.cc
  newmarkFracDamp.cc
  acouCombustion.cc
  smoothPDE.cc
  SinglePDE.cc
  newmark.cc
  basePDE.cc
  bdf2.cc
  trapezoidal.cc
  pseudoTS.cc
  timestepping.cc
  mixedEqnMap.cc
  extLBMPDE.cc
  )
# mpcciPDE.cc

ADD_LIBRARY(pde STATIC ${PDE_SRCS})

SET(TARGET_LL
  algsys-olas
  cfsgeneral
  coupledpde
  datainout
  design
  domain
  driver
  elements
  forms
  logging
  materials
  matvec
  optimization
  paramh
  scripting
  utils
)

IF(USE_HDF5)
  SET(TARGET_LL
    ${TARGET_LL}
    siminouthdf5
    )
ENDIF(USE_HDF5)

IF(USE_MPCCI)
  SET(TARGET_LL
    ${TARGET_LL}
    mpccicpl
    ${MPCCI_LIBRARIES}
    )
ENDIF(USE_MPCCI)

TARGET_LINK_LIBRARIES(pde ${TARGET_LL})