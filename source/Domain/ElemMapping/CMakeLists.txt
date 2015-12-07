SET(MAPPINGS_SRCS 
  EdgeFace.cc
  Elem.cc
  ElemShapeMap.cc
  EntityLists.cc
  )


ADD_LIBRARY(elemmapping STATIC ${MAPPINGS_SRCS})

SET(TARGET_LL
  utils
  domain
  )


TARGET_LINK_LIBRARIES(elemmapping
  ${TARGET_LL}
  )
  
