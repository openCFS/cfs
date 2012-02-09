# ------------- MAKEFILE FOR FILEREADER FOR COUPLING -----------------------
# Author: LSE.
# Date:   05.09.2003.
# --------------------------------------------------------------------------
#

SET(CPLREADER_ANSYS_SRCS
  FileReader_ANSYS.cc
  FileReader_MKHDF5.cc
  FileReader_NACS.cc
)

ADD_LIBRARY(cplreader_ansys STATIC ${CPLREADER_ANSYS_SRCS})

TARGET_LINK_LIBRARIES(cplreader_ansys
  cfsgeneral
  ${BOOST_REGEX_LIB}
  ${BOOST_FILESYSTEM_LIB}
)