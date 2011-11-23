SET(ODESOLVER_SRCS ODESolver_RKF45.cc ODESolver_ExplEuler.cc 
  ODESolver_Rosenbrock.cc piezoSwitch.cc )

ADD_LIBRARY(odesolver STATIC ${ODESOLVER_SRCS})
