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
  AcousticMixedPDE.cc
  ElecCurrentPDE.cc
  TestPDE.cc
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
