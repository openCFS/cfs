IF(USE_XERCES)
  INCLUDE_DIRECTORIES(${XERCES_INCLUDE_DIR})
ENDIF(USE_XERCES)

IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
ENDIF(USE_IPOPT)

IF(USE_XERCES)
  INCLUDE_DIRECTORIES(${XERCES_INCLUDE_DIR})
ENDIF(USE_XERCES)


SET(DOMAIN_SRCS 
  #ansatzFct.cc
  bcs.cc
  domain.cc
  edgeFace.cc
  elem.cc
  entityList.cc
  grid.cc
  resultInfo.cc
  )

ADD_LIBRARY(domain STATIC ${DOMAIN_SRCS})

SET(TARGET_LL
  elements
#s  optimization
  cfsgeneral
  driver
  utils
  pde
  coupledpde
  gridcfs
  scripting
  matvec
# DIRTY, DIRTYYYY HACK!
  datainout
  simouttext
  datainout
)

IF(USE_INTERPOLATION)
  SET(TARGET_LL
    ${TARGET_LL}
    ${CGAL_LIBRARY}
    ${BOOST_THREAD_LIB}
  )
ENDIF(USE_INTERPOLATION)

IF(USE_XERCES)
  SET(TARGET_LL
    ${TARGET_LL}
    ${XERCES_LIBRARY}
  )
ENDIF(USE_XERCES)


TARGET_LINK_LIBRARIES(domain
  ${TARGET_LL}
  )
