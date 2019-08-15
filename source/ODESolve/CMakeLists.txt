SET(ODESOLVER_SRCS 
ODESolver_RKF45.cc 
ODESolver_ExplEuler.cc 
ODESolver_Rosenbrock.cc
PiezoSwitch.cc)

ADD_LIBRARY(odesolver STATIC ${ODESOLVER_SRCS})

ADD_DEPENDENCIES(odesolver boost)
