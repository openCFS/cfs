IF(USE_MPCCI)
  INCLUDE_DIRECTORIES(${MPCCI_INCLUDE_DIR})
ENDIF(USE_MPCCI)

IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

SET(PDE_SRCS
  acouCombustion.cc
  acousticMixedPDE.cc
  acousticPDE.cc
  basePDE.cc
  bdf2.cc
  disContEqnMap.cc
  elecPDE.cc
  eqnMap.cc
  fluidMechPDE.cc
  heatCondPDE.cc
  magEdgePDE.cc
  magneticPDE.cc
  mechPDE.cc
  mixedEqnMap.cc
  mpcciPDE.cc
  newmark.cc
  newmarkFracDamp.cc
  newmarkFracDampMech.cc
  pdememento.cc
  pseudoTS.cc
  SinglePDE.cc
  smoothPDE.cc
  StdPDE.cc
  timestepping.cc
  trapezoidal.cc
  )

ADD_LIBRARY(pde STATIC ${PDE_SRCS})

SET(TARGET_LL
  datainout
  domain
  algsys-olas
  forms
  driver
  paramh
  matvec
  utils
  ${BOOST_SERIALIZATION_HDF5_LIB}
)

IF(USE_MPCCI)
  SET(TARGET_LL
    ${TARGET_LL}
    mpccicpl
    ${MPCCI_LIBRARIES}
    )
ENDIF(USE_MPCCI)

TARGET_LINK_LIBRARIES(pde ${TARGET_LL})
