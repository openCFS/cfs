SET(SIMOUTUNV_SRCS
  simOutUnv.cc
  simInputUnv.cc
  unv_reader.cc
  unv.cc
  unv_dat.cc
  unv_if.cc
)

ADD_LIBRARY(simoutunv STATIC ${SIMOUTUNV_SRCS})
