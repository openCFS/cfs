SET(UTILS_SRCS
  math/givensrotation.cc
  math/iterativerefinement.cc
)


IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(UTILS_SRCS
    ${UTILS_SRCS}
    math/croutlu.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(utils-olas STATIC ${UTILS_SRCS})

TARGET_LINK_LIBRARIES(utils-olas
  cfsgeneral
  lapack-olas
  matvec
  precond-olas
  solver-olas
  utils-olas
)
