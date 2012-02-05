SET(PARAMH_SRCS
  ParamNode.cc
  XMLSkeletonConf.cc
  CFSOLASParams.cc
)

SET(TARGET_LL
  ${BOOST_IOSTREAMS_LIB}
  ${BOOST_FILESYSTEM_LIB}
  ${BOOST_SYSTEM_LIB}
  cfsgeneral
  datainout
  domain
  driver
  elements
  graph-olas
  matvec
  pde
  precond-olas
  solver-olas
  utils
)

IF(USE_XERCES)

  INCLUDE_DIRECTORIES(${XERCES_INCLUDE_DIR})

  SET(PARAMH_SRCS
    ${PARAMH_SRCS}
    Xerces.cc)

  SET(TARGET_LL ${TARGET_LL} ${XERCES_LIBRARY})

ENDIF(USE_XERCES)

ADD_LIBRARY(paramh STATIC ${PARAMH_SRCS})

TARGET_LINK_LIBRARIES(paramh utils ${TARGET_LL})
