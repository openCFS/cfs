SET(SIMOUTINFO_SRCS SimOutputInfo.cc)

ADD_LIBRARY(simoutinfo STATIC ${SIMOUTINFO_SRCS})

TARGET_LINK_LIBRARIES(simoutinfo 
  datainout
  domain
  matvec
  paramh
  pde
  utils
)