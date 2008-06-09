SET(UTILS_SRCS
  environment-olas.cc
  output.cc
  tools-olas.cc
  profiler-olas.cc
  math/givensrotation.cc
  math/iterativerefinement.cc)

ADD_LIBRARY(utils-olas STATIC ${UTILS_SRCS})

TARGET_LINK_LIBRARIES(utils-olas
  lapack-olas)
