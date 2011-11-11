SET(UTILS_SRCS
    math/GivensRotation.cc
    math/IterativeRefinement.cc)


IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(UTILS_SRCS
    ${UTILS_SRCS}
    math/CroutLU.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(utils-olas STATIC ${UTILS_SRCS})

TARGET_LINK_LIBRARIES(utils-olas
  matvec
)
