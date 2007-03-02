# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(ODESOLVER_SRCS ODESolver_RKF45.cc ODESolver_ExplEuler.cc 
  ODESolver_Rosenbrock.cc)

ADD_LIBRARY(odesolver STATIC ${ODESOLVER_SRCS})
