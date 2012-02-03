SET(SIMOUTUNV_SRCS
  SimOutputUnv.cc
  SimInputUnv.cc
  unv.cc
  unv_dat.cc
  unv_if.cc
  unv_reader.cc
)

ADD_LIBRARY(simoutunv STATIC ${SIMOUTUNV_SRCS})

TARGET_LINK_LIBRARIES(simoutunv
  cfsgeneral
  datainout
  domain
  elements
  logging
  matvec
  paramh
  pde
  utils
)