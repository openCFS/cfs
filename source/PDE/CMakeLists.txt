# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

IF(USE_MPCCI)
  INCLUDE_DIRECTORIES(${MPCCI_INCLUDE_DIR})
ENDIF(USE_MPCCI)

IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

SET(PDE_SRCS basePDE.cc StdPDE.cc SinglePDE.cc elecPDE.cc mechPDE.cc
  acousticPDE.cc newmark.cc timestepping.cc smoothPDE.cc magneticPDE.cc
  trapezoidal.cc pdememento.cc acouFlowNoise.cc newmarkFracDamp.cc 
  newmarkFracDampMech.cc mpcciPDE.cc heatCondPDE.cc 
  acousticXYZPDE.cc stokesFluidPDE.cc eqnMap.cc bubblePDE.cc 
  acouCombustion.cc acousticMixedPDE.cc magEdgePDE.cc)

ADD_LIBRARY(pde STATIC ${PDE_SRCS})

SET(TARGET_LL
  datainout
  domain
  algsys-olas
  forms
  driver
  utils
  odedescr
  odesolver
  paramh
)

IF(USE_MPCCI)
  SET(TARGET_LL
    ${TARGET_LL}
    mpccicpl
    ${MPCCI_LIBRARIES}
    )
ENDIF(USE_MPCCI)

TARGET_LINK_LIBRARIES(pde ${TARGET_LL})
