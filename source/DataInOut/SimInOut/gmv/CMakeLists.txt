IF(USE_GMV_INPUT)

SET(SIMINGMV_SRCS
  gmvread.cc
  simInputGMV.cc
)

  ADD_LIBRARY(simingmv STATIC ${SIMINGMV_SRCS})

  TARGET_LINK_LIBRARIES(simingmv
    cfsgeneral
    domain
    logging
    matvec
    paramh
    pde
    utils
  )

ENDIF(USE_GMV_INPUT)

IF(USE_GMV_OUTPUT)

  SET(SIMOUTGMV_SRCS simOutGMV.cc)

  ADD_LIBRARY(simoutgmv STATIC ${SIMOUTGMV_SRCS})

  TARGET_LINK_LIBRARIES(simoutgmv
    cfsgeneral
    datainout
    domain
    logging
    matvec
    paramh
    pde
    utils
  )

ENDIF(USE_GMV_OUTPUT)
