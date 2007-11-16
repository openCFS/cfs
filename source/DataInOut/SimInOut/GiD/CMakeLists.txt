SET(SIMOUTGID_SRCS
  simOutGiD.cc
  gidpost/gidpost.cpp
#  gidpost/gidpostfor.c
  gidpost/gidpostInt.cpp
  gidpost/zlibint1.c
  gidpost/zlibint2.c)

ADD_DEFINITIONS("-Df2cFortran")

ADD_LIBRARY(simoutgid STATIC ${SIMOUTGID_SRCS})

TARGET_LINK_LIBRARIES(simoutgid)
