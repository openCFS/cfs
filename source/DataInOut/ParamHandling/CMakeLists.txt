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

# there is a problem in ParamNode.cc and Xerces.cc and gcc 4.5.1 with 
# boost/iostreams/filter/symmetric.hpp:126:13: error: suggest parentheses around ‘&&’ within ‘||’
# boost/iostreams/detail/streambuf/indirect_streambuf.hpp:285:5: error: suggest parentheses around ‘&&’ within ‘||’
# TODO FIXME check if this is still necessary on FE-Space with new boost
SET_SOURCE_FILES_PROPERTIES(${PARAMH_SRCS}
  PROPERTIES COMPILE_FLAGS "-Wno-parentheses")

ADD_LIBRARY(paramh STATIC ${PARAMH_SRCS})

TARGET_LINK_LIBRARIES(paramh utils ${TARGET_LL})
