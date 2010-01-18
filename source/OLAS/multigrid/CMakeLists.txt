SET(MATVEC_SRCS stdmatrix.cc stdvector.cc sbmvector.cc sbmmatrix.cc
	        generatematvec.cc)

ADD_DEFINITIONS()

ADD_LIBRARY(matvec STATIC ${MATVEC_SRCS})


