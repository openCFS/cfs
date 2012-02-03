SET(ODESOLVER_SRCS
  ODESolver_Rosenbrock.cc
  ODESolver_RKF45.cc
  ODESolver_ExplEuler.cc
  piezoSwitch.cc
 )

ADD_LIBRARY(odesolver STATIC ${ODESOLVER_SRCS})

TARGET_LINK_LIBRARIES(odesolver
  cfsgeneral
  matvec
  utils
)
