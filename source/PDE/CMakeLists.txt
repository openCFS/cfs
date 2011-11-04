IF(USE_MPCCI)
  INCLUDE_DIRECTORIES(${MPCCI_INCLUDE_DIR})
ENDIF(USE_MPCCI)

IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

SET(PDE_SRCS
#  acouCombustion.cc
# acousticMixedPDE.cc
  basePDE.cc
#  disContEqnMap.cc
  elecPDE.cc
#  eqnMap.cc
#  fluidMechPDE.cc
#  heatCondPDE.cc
  magEdgePDE.cc
#  magneticPDE.cc
# magneticScalarPDE.cc
#  mechPDE.cc
#  mixedEqnMap.cc
  newmark.cc
  pseudoTS.cc
  SinglePDE.cc
#  smoothPDE.cc
  StdPDE.cc
  timestepping.cc
  trapezoidal.cc
  )

ADD_LIBRARY(pde STATIC ${PDE_SRCS})

SET(TARGET_LL
  datainout
  domain
  algsys-olas
  driver
  mathparser
  paramh
  matvec
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
    )
ENDIF(USE_MPCCI)

TARGET_LINK_LIBRARIES(pde ${TARGET_LL})
