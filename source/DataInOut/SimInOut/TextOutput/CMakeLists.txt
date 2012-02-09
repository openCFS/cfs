SET(SIMOUTTEXT_SRCS textSimOutput.cc)

ADD_LIBRARY(simouttext STATIC ${SIMOUTTEXT_SRCS})

TARGET_LINK_LIBRARIES(simouttext
  cfsgeneral
  datainout
  domain
  elements
  matvec
  paramh
  pde
  utils
)
