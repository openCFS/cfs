SET(GENERAL_SRCS
  environment.cc
  exception.cc
)

ADD_LIBRARY(cfsgeneral STATIC ${GENERAL_SRCS})

TARGET_LINK_LIBRARIES(cfsgeneral
  logging
  matvec
)

