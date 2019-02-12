SET(SIMINREFELEMS_SRCS SimInputRefElems.cc)

ADD_LIBRARY(siminrefelems STATIC ${SIMINREFELEMS_SRCS})

ADD_DEPENDENCIES(siminrefelems boost)

IF(TARGET cgal)
  ADD_DEPENDENCIES(siminrefelems cgal)
ENDIF()
