SET(CPLREADER_ENSIGHT_SRCS
  FileReader_EnSight.cc
  )

ADD_LIBRARY(cplreader_ensight STATIC ${CPLREADER_ENSIGHT_SRCS})

TARGET_LINK_LIBRARIES(cplreader_ensight
  cplreader_vtkmb
  domain
  ${VTK_LIB}
)

SET_SOURCE_FILES_PROPERTIES(${CPLREADER_ENSIGHT_SRCS} MergePointExample.cc
  PROPERTIES COMPILE_FLAGS "-Wno-deprecated")

ADD_EXECUTABLE(MergePointExample MergePointExample.cc)

TARGET_LINK_LIBRARIES(MergePointExample ${VTK_LIB} -pthread)
