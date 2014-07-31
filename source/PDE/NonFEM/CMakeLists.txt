SET(NON_FEM_SRCS
  LatticeBoltzmann.cc
  )

ADD_LIBRARY(non_fem_pde STATIC ${NON_FEM_SRCS})

SET(TARGET_LL
  utils
)


TARGET_LINK_LIBRARIES(non_fem_pde ${TARGET_LL})