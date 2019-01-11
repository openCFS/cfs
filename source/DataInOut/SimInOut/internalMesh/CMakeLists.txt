SET(INTERNALMESH_SRCS InternalMesh.cc)

ADD_LIBRARY(internalmesh STATIC ${INTERNALMESH_SRCS})

ADD_DEPENDENCIES(internalmesh boost)

IF(TARGET cgal)
  ADD_DEPENDENCIES(internalmesh cgal)
ENDIF()
