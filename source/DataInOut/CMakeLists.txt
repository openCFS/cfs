SET(DATAINOUT_SRCS 
  freqfunc.cc
  XMLMaterialHandler.cc
  PlainMaterialHandler.cc
  WriteInfo.cc
  resultHandler.cc
  simOutput.cc
  postProc.cc)

ADD_LIBRARY(datainout STATIC ${DATAINOUT_SRCS})
