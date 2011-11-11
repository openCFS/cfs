SET(MAPPINGS_SRCS 
  EdgeFace.cc
  Elem.cc
  ElemShapeMap.cc
  EntityLists.cc
  Jacobian.cc
  )


ADD_LIBRARY(elemmapping STATIC ${MAPPINGS_SRCS})

SET(TARGET_LL
  matvec
  cfsgeneral
  utils
  mesh
  febasis
  )


TARGET_LINK_LIBRARIES(elemmapping
  ${TARGET_LL}
  )
  