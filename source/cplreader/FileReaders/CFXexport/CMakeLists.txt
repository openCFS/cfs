SET(CPLREADER_CFXEXPORT_SRCS
  FileReader_CFXexport.cc
)

ADD_LIBRARY(cplreader_cfxexport STATIC ${CPLREADER_CFXEXPORT_SRCS})

TARGET_LINK_LIBRARIES(cplreader_cfxexport
  cfsgeneral
)