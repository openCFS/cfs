SET(DATAINOUT_SRCS 
  freqfunc.cc
  XMLMaterialHandler.cc
  WriteInfo.cc
  resultHandler.cc
  simOutput.cc
  postProc.cc)

ADD_LIBRARY(datainout STATIC ${DATAINOUT_SRCS})

SET(TARGET_LL
  materials
  simouttext  
)

TARGET_LINK_LIBRARIES(datainout ${TARGET_LL})
