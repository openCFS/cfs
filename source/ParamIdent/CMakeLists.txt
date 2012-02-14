SET(PARAMIDENT_SRCS
  piezoParamIdent.cc
  piezoParamIdent_methods.cc
  nonlinLandweber.cc
  nuMethods.cc
  leastSquare.cc
)

ADD_LIBRARY(paramident STATIC ${PARAMIDENT_SRCS})

TARGET_LINK_LIBRARIES(paramident
  cfsgeneral
  coupledpde
  datainout
  domain
  driver
  lapack-olas
  materials
  matvec
  paramh
  pde
  utils
)