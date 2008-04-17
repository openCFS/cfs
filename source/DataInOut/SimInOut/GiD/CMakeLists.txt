SET(SIMOUTGID_SRCS
  simOutGiD.cc
  gidpost/gidpost.cpp
#  gidpost/gidpostfor.c
  gidpost/gidpostInt.cpp)

ADD_DEFINITIONS("-Df2cFortran")

ADD_LIBRARY(simoutgid STATIC ${SIMOUTGID_SRCS})