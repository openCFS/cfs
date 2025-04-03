set(GINKGO_SRCS GinkgoSolver.cc)

add_library(ginkgo-olas STATIC ${GINKGO_SRCS})

target_link_libraries(ginkgo-olas ${GINKGO_LIBRARY})

if(USE_CUDA)
  # this works manually: tailing -L"/usr/local/cuda-12.8/targets/x86_64-linux/lib/stubs" -L"/usr/local/cuda-12.8/targets/x86_64-linux/lib" -lcudart -lcusparse -lcublas -lcurand 
  #cmake_print_variables(CMAKE_CUDA_HOST_IMPLICIT_LINK_DIRECTORIES) # /usr/local/cuda-12.8/targets/x86_64-linux/lib/stubs;/usr/local/cuda-12.8/targets/x86_64-linux/lib
  assert_set(CMAKE_CUDA_HOST_IMPLICIT_LINK_DIRECTORIES)
  target_link_directories(ginkgo-olas PUBLIC ${CMAKE_CUDA_HOST_IMPLICIT_LINK_DIRECTORIES})
  target_link_libraries(ginkgo-olas "$<LINK_GROUP:RESCAN,cudart;cusparse;cublas;curand>")
endif()

add_dependencies(ginkgo-olas cfsdeps)
