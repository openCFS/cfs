SET(ALGSYS_SRCS
    algebraicSys.cc
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

SET(TARGET_LL
  matvec
  solver-olas
  precond-olas
  graph-olas
)

TARGET_LINK_LIBRARIES(algsys-olas ${TARGET_LL})
