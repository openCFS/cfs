IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

SET(DOMAIN_SRCS domain.cc grid.cc entityList.cc resultInfo.cc
  ansatzFct.cc bcs.cc edgeFace.cc)

ADD_LIBRARY(domain STATIC ${DOMAIN_SRCS})
