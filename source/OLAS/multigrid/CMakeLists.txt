SET(MATVEC_SRCS StdMatrix.cc StdVector.cc SBM_Vector.cc SBM_Matrix.cc
	        generatematvec.cc)

ADD_DEFINITIONS()

ADD_LIBRARY(matvec STATIC ${MATVEC_SRCS})


