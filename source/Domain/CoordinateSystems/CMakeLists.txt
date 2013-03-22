SET(COORDSYS_SRCS
  CartesianCoordSystem.cc
  CoordSystem.cc
  CylCoordSystem.cc
  DefaultCoordSystem.cc
  PolarCoordSystem.cc
  TrivialCartesianCoordSystem.cc
  )

ADD_LIBRARY(coordsystems STATIC ${COORDSYS_SRCS})

SET(TARGET_LL
  matvec
  utils
  paramh
  mesh
  cfsgeneral
  )

TARGET_LINK_LIBRARIES(coordsystems
  ${TARGET_LL}
  )
