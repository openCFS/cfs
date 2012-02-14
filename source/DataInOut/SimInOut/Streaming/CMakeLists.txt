SET(SIMOUTSTREAMING_SRCS SimOutputStreaming.cc)

ADD_LIBRARY(simoutstreaming STATIC ${SIMOUTSTREAMING_SRCS})

TARGET_LINK_LIBRARIES(simoutstreaming
  cfsgeneral
  datainout
  domain
  matvec
  paramh
  pde
  utils
)
