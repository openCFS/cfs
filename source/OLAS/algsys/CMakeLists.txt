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
