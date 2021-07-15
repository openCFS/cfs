SET(LBMSOLVER_SRCS
  LatticeBoltzmann.cc
  )

ADD_LIBRARY(lbm_solver STATIC ${LBMSOLVER_SRCS})

SET(TARGET_LL
  utils optimization
)

TARGET_LINK_LIBRARIES(lbm_solver ${TARGET_LL})
