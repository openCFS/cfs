# ------------- MAKEFILE FOR FILEREADER FOR COUPLING -----------------------
# Author: LSE.
# Date:   05.09.2003.
# --------------------------------------------------------------------------
#

IF(CPLREADER_REBUILD_DOCU)
  SUBDIRS(latex)
  SET(PDF_FILE "${CMAKE_CURRENT_BINARY_DIR}/latex/cplreader_docu.pdf")
  SET(WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/latex")
ELSE(CPLREADER_REBUILD_DOCU)
  SET(PDF_FILE "${CMAKE_CURRENT_SOURCE_DIR}/latex/cplreader_docu.pdf")
  SET(WORK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/latex")
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
  ADD_DEPENDENCIES(cplreader_docu cplreader_docu_pdf)

  ADD_CUSTOM_TARGET(cplreader_docu_copy ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CFS_BINARY_DIR}/source/cplreader/Documentation/latex/cplreader_docu.pdf ${CFS_SOURCE_DIR}/source/cplreader/Documentation/latex/cplreader_docu.pdf
  COMMENT "Copying cplreader_docu.pdf back to ${CFS_SOURCE_DIR}/source/cplreader/Documentation..."
  )

  ADD_DEPENDENCIES(cplreader_docu_copy cplreader_docu_pdf)
ENDIF(CPLREADER_REBUILD_DOCU)
