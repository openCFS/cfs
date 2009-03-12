#-------------------------------------------------------------------------------
# CFX specific settings
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES(${CFX_V110_INCLUDE_DIR})

SET(CPLREADER_CFX_SRCS
  FileReader_CFX.cc
  glibchack.cc
  OpenFile.F
  CloseFile.F
  iocntf_wrap.F
  RedHdr.F
  ReadLong.F
  RedSht.F
  tools.F
)

ADD_LIBRARY(cplreader_cfx STATIC ${CPLREADER_CFX_SRCS})

SET_TARGET_PROPERTIES(cplreader_cfx
  PROPERTIES
  LINKER_LANGUAGE CXX
)

TARGET_LINK_LIBRARIES(cplreader_cfx
  ${CFX_V110_LIBIO} # CFX I/O library
  ${BOOST_REGEX_LIB}
)
