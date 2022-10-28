SET(UTILS_SRCS
    math/GivensRotation.cc
    math/IterativeRefinement.cc)

SET(UTILS_SRCS
    ${UTILS_SRCS}
    math/CroutLU.cc
    )

ADD_LIBRARY(utils-olas STATIC ${UTILS_SRCS})

TARGET_LINK_LIBRARIES(utils-olas
  matvec
)
