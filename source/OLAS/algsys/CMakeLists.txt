SET(ALGSYS_SRCS
  standardsys.cc
  sbmsystem.cc
  basesystem.cc
  baseentrymanipulator.cc
  generateidbchandler.cc
    )

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(ALGSYS_SRCS
    ${ALGSYS_SRCS}
    idbchandler.cc
    idbchandlerpenalty.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(algsys-olas STATIC ${ALGSYS_SRCS})

TARGET_LINK_LIBRARIES(algsys-olas
  cfsgeneral
  datainout
  graph-olas
  logging
  matvec
  paramh
  precond-olas
  solver-olas
  utils
)

