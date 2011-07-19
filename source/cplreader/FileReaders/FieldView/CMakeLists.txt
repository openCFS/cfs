SET(CPLREADER_FV_SRCS
  FileReader_fvuns.cc
  readfvuns.cpp
  )

ADD_LIBRARY(cplreader_fieldview STATIC ${CPLREADER_FV_SRCS})

#ADD_EXECUTABLE(readfvuns readfvuns.cpp)
