SET(CPLREADER_FV_SRCS
  FileReader_fvuns.cc
  FieldViewFile.cpp
  )

ADD_LIBRARY(cplreader_fieldview STATIC ${CPLREADER_FV_SRCS})

TARGET_LINK_LIBRARIES(cplreader_fieldview
  domain
  matvec
  ${BOOST_FILESYSTEM_LIB}
)
