SET(SIMINMESH_SRCS simInputMESH.cc)

ADD_LIBRARY(siminmesh STATIC ${SIMINMESH_SRCS})

TARGET_LINK_LIBRARIES(siminmesh
  cfsgeneral
  domain
  matvec
  paramh
  utils
)