SET(MATVEC_SRCS
  generatematvec.cc
  stdmatrix.cc
  sbmmatrix.cc
  sbmvector.cc
  stdvector.cc
)

ADD_DEFINITIONS()

ADD_LIBRARY(matvec STATIC ${MATVEC_SRCS})


