# Make sure the compiler can find include files.
# INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})

SET(SIMOUTENSIGHT_SRCS
  simOutputEnsightGold.cc
)

ADD_LIBRARY(simoutensight STATIC ${SIMOUTENSIGHT_SRCS})

