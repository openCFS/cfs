# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(PARAMH_SRCS PlainXMLParamHandler.cc XMLSkeletonConf.cc CFSOLASParams.cc
  ParamNode.cc)

IF(USE_XERCES)
  SET(PARAMH_SRCS ${PARAMH_SRCS} XMLParamHandler.cc XMLParserErrorHandler.cc
    Xerces.cc)
ENDIF(USE_XERCES)

ADD_LIBRARY(paramh STATIC ${PARAMH_SRCS})

IF(USE_XERCES)
  TARGET_LINK_LIBRARIES(paramh xerces-c)
ENDIF(USE_XERCES)
