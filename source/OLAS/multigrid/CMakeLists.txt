IF(EXPLICIT_TEMPLATE_INSTANTIATION)

SET(MG_SRCS
    ${MG_SRCS}
    amg.cc
    hierarchylevel.cc
    topology.cc
    transfer.cc
    gaussseidel.cc
    jacobi.cc
    smoother.cc
    transfer.cc
    depgraph.cc
    prematrix.cc
    agglomerate.cc
    AFWsmoother.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)


ADD_LIBRARY(multigrid-olas STATIC ${MG_SRCS})

SET(TARGET_LL
  multigrid-olas
)

TARGET_LINK_LIBRARIES(multigrid-olas ${TARGET_LL})