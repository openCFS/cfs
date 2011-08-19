SET(SCRIPTING_SRCS scriptable.cc cfsmessenger.cc)

ADD_LIBRARY(scripting STATIC ${SCRIPTING_SRCS})
