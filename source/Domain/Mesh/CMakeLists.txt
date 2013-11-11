IF(USE_XERCES)
  INCLUDE_DIRECTORIES(${XERCES_INCLUDE_DIR})
ENDIF(USE_XERCES)


SET(MESH_SRCS 
  Grid.cc
  NcInterfaces/BaseNcInterface.cc
  NcInterfaces/MortarInterface.cc
  GridCFS/GridCFS.cc
  )


ADD_LIBRARY(mesh STATIC ${MESH_SRCS})

SET(TARGET_LL
  utils
  elemmapping
  coordsystems
  datainout
  paramh
  matvec
  cfsgeneral
  integration
  )

IF(USE_CGAL)
  SET(TARGET_LL
    ${TARGET_LL}
    ${CGAL_LIBRARY}
    ${BOOST_THREAD_LIB}
  )
ENDIF(USE_CGAL)

IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC" AND USE_OPENMP)
  SET(TARGET_LL
    ${TARGET_LL}
    ${BOOST_THREAD_LIB}
  )
ENDIF()

IF(USE_XERCES)
  SET(TARGET_LL
    ${TARGET_LL}
    ${XERCES_LIBRARY}
  )
ENDIF(USE_XERCES)

LIST(REMOVE_DUPLICATES TARGET_LL)

TARGET_LINK_LIBRARIES(mesh
  ${TARGET_LL}
  )
