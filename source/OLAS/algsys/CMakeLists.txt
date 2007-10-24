# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(ALGSYS_SRCS basesystem.cc standardsys.cc sbmsystem.cc olascomm.cc olasparams.cc
      generateidbchandler.cc generateentrymanipulator.cc)

ADD_DEFINITIONS()

ADD_LIBRARY(algsys-olas STATIC ${ALGSYS_SRCS})

SET(TARGET_LL
  matvec-olas
  solver-olas
  precond-olas
)

TARGET_LINK_LIBRARIES(algsys-olas ${TARGET_LL})
