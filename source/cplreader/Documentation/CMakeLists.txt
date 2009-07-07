# ------------- MAKEFILE FOR FILEREADER FOR COUPLING -----------------------
# Author: LSE.
# Date:   05.09.2003.
# --------------------------------------------------------------------------
#

# If documentation is to be rebuilt set this to 1
SET(BUILD_CPLREADER_DOCU 0)

IF(BUILD_CPLREADER_DOCU)
  SUBDIRS(latex)
  SET(PDF_FILE "${CMAKE_CURRENT_BINARY_DIR}/latex/cplreader_docu.pdf")
  SET(WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/latex")
ELSE(BUILD_CPLREADER_DOCU)
  SET(PDF_FILE "${CMAKE_CURRENT_SOURCE_DIR}/latex/cplreader_docu.pdf")
  SET(WORK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/latex")
ENDIF(BUILD_CPLREADER_DOCU)

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

IF(BUILD_CPLREADER_DOCU)
  ADD_DEPENDENCIES(cplreader_docu cplreader_docu_pdf)
ENDIF(BUILD_CPLREADER_DOCU)
