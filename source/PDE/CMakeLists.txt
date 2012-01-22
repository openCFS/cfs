IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

SET(PDE_SRCS
  BasePDE.cc
  ElecPDE.cc
  MagEdgePDE.cc
  SinglePDE.cc
  StdPDE.cc
  HeatPDE.cc
  AcousticPDE.cc
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
)

IF(USE_HDF5)
  SET(TARGET_LL
    ${TARGET_LL}
    siminouthdf5
    )
ENDIF(USE_HDF5)

TARGET_LINK_LIBRARIES(pde ${TARGET_LL})
