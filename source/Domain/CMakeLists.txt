IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

SET(DOMAIN_SRCS 
  ansatzFct.cc
  bcs.cc
  domain.cc
  edgeFace.cc
  entityList.cc
  grid.cc
  resultInfo.cc
  )

ADD_LIBRARY(domain STATIC ${DOMAIN_SRCS})

SET(TARGET_LL
  elements
  optimization
  general
  driver
  utils
  pde
  coupledpde
  gridcfs
  scripting
# DIRTY, DIRTYYYY HACK!
  datainout
  simouttext
  datainout
)

IF(USE_INTERPOLATION)
  SET(TARGET_LL
    ${TARGET_LL}
    CGAL
  )
ENDIF(USE_INTERPOLATION)

TARGET_LINK_LIBRARIES(domain
  ${TARGET_LL}
  )