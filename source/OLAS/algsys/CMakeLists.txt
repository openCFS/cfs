# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(ALGSYS_SRCS
    basesystem.cc
    baseentrymanipulator.cc
    generateidbchandler.cc
    olascomm.cc
    olasparams.cc
    sbmsystem.cc
    standardsys.cc
    )

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(ALGSYS_SRCS
    ${ALGSYS_SRCS}
    idbchandler.cc
    idbchandlerpenalty.cc
    olascomm_templ_impl.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(algsys-olas STATIC ${ALGSYS_SRCS})

SET(TARGET_LL
  matvec
  solver-olas
  precond-olas
  graph-olas
)

TARGET_LINK_LIBRARIES(algsys-olas ${TARGET_LL})
