SET(PARAMIDENT_SRCS piezoParamIdent.cc piezoParamIdent_methods.cc nuMethods.cc
  nonlinLandweber.cc 
  leastSquare.cc)

ADD_LIBRARY(paramident STATIC ${PARAMIDENT_SRCS})
