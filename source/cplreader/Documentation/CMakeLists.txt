# ------------- MAKEFILE FOR FILEREADER FOR COUPLING -----------------------
# Author: LSE.
# Date:   05.09.2003.
# --------------------------------------------------------------------------
#

SET(DOC_DEPS index.html logo.png)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/GetHTMLDocumentation.cc
  COMMAND perl ${CMAKE_CURRENT_SOURCE_DIR}/docfilestocpp.pl ${CMAKE_CURRENT_BINARY_DIR}/GetHTMLDocumentation.cc GetHTMLDocumentation '*.html *.png'
  MAIN_DEPENDENCY docfilestocpp.pl
  DEPENDENCIES ${DOC_DEPS}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Generating C++ source from cplreader docs...")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CPLREADER_DOC_SRCS 
  ${CMAKE_CURRENT_BINARY_DIR}/GetHTMLDocumentation.cc
  WriteDocumentation.cc
)

ADD_LIBRARY(cplreader_doc STATIC ${CPLREADER_DOC_SRCS})

