# ------------- MAKEFILE FOR FILEREADER FOR COUPLING -----------------------
# Author: LSE.
# Date:   05.09.2003.
# --------------------------------------------------------------------------
#

IF(CPLREADER_REBUILD_DOCU)
  SET(PDF_FILE "${CFS_BINARY_DIR}/share/doc/user/cplreader/cplreader_docu.pdf")
ELSE(CPLREADER_REBUILD_DOCU)
  SET(PDF_FILE "${CFS_SOURCE_DIR}/share/doc/user/cplreader/cplreader_docu.pdf")
ENDIF(CPLREADER_REBUILD_DOCU)


ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/GetPDFDocumentation.cc
                   COMMAND perl ${CMAKE_CURRENT_SOURCE_DIR}/docfilestocpp.pl ${CMAKE_CURRENT_BINARY_DIR}/GetPDFDocumentation.cc GetPDFDocumentation '*.pdf'
                   DEPENDS docfilestocpp.pl ${PDF_FILE}
                   WORKING_DIRECTORY ${WORK_DIR}
                   COMMENT "Generating C++ source from cplreader docs...")
                     
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CPLREADER_DOC_SRCS 
  ${CMAKE_CURRENT_BINARY_DIR}/GetPDFDocumentation.cc
  WriteDocumentation.cc
)

ADD_LIBRARY(cplreader_docu STATIC ${CPLREADER_DOC_SRCS})

IF(CPLREADER_REBUILD_DOCU)
  ADD_DEPENDENCIES(cplreader_docu doc-user)
ENDIF(CPLREADER_REBUILD_DOCU)
