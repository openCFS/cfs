IF(USE_XERCES)
  INCLUDE_DIRECTORIES(${XERCES_INCLUDE_DIR})
ENDIF(USE_XERCES)


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

IF(USE_XERCES)
  SET(TARGET_LL
    ${TARGET_LL}
    ${XERCES_LIBRARY}
    )
ENDIF(USE_XERCES)

TARGET_LINK_LIBRARIES(datainout ${TARGET_LL})
