SET(SIMOUTUNV_SRCS
  SimOutputUnv.cc
  SimInputUnv.cc
  unv_reader.cc
  unv.cc
  unv_dat.cc
  unv_if.cc
)

ADD_LIBRARY(simoutunv STATIC ${SIMOUTUNV_SRCS})

ADD_DEPENDENCIES(simoutunv boost)

IF(TARGET cgal)
  ADD_DEPENDENCIES(simoutunv cgal)
ENDIF()

