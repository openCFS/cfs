SET(UTILS_SRCS
    math/givensrotation.cc
    math/iterativerefinement.cc)

IF(PROFILING)
  SET(UTILS_SRCS
    ${UTILS_SRCS}
    profiler-olas.cc
    )
ENDIF(PROFILING)

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(UTILS_SRCS
    ${UTILS_SRCS}
    math/croutlu.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(utils-olas STATIC ${UTILS_SRCS})

TARGET_LINK_LIBRARIES(utils-olas
  matvec
)
