# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(PARAMIDENT_SRCS piezoParamIdent.cc piezoParamIdent_methods.cc nuMethods.cc
  nonlinLandweber.cc 
  leastSquare.cc)

ADD_LIBRARY(paramident STATIC ${PARAMIDENT_SRCS})
