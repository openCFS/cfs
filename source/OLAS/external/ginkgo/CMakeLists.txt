set(GINKGO_SRCS GinkgoSolver.cc)

add_library(ginkgo-olas STATIC ${GINKGO_SRCS})

target_link_libraries(ginkgo-olas ${GINKGO_LIBRARY})

add_dependencies(ginkgo-olas cfsdeps)
