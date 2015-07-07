IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

INCLUDE_DIRECTORIES(${XERCES_INCLUDE_DIR})

SET(DOMAIN_SRCS 
  BCs.cc
  Domain.cc
  )

ADD_LIBRARY(domain STATIC ${DOMAIN_SRCS})

SET(TARGET_LL
  cfsgeneral
  optimization
  mesh
  coordsystems
  results
  pde
  coupledpde
  driver
  mathparser
  utils
# DIRTY, DIRTYYYY HACK! WHY?????
  datainout
  simouttext
  datainout
)

SET(TARGET_LL
  ${TARGET_LL}
  ${XERCES_LIBRARY}
)

TARGET_LINK_LIBRARIES(domain
  ${TARGET_LL}
  )
