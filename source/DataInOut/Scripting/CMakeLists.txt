SET(SCRIPTING_SRCS
  scriptable.cc
  cfsmessenger.cc
)

ADD_LIBRARY(scripting STATIC ${SCRIPTING_SRCS})

TARGET_LINK_LIBRARIES(scripting
  cfsgeneral
  datainout
  domain
  driver
  logging
  matvec
  pde
  utils
)
