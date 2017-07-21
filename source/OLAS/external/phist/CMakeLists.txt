SET(PHIST_SRCS
  PhistEigenSolver.cc
)

include_directories(/usr/lib64/mpi/gcc/openmpi/include
                    /home/fwein/code/ghost/build/hwloc-inst/include)

ADD_LIBRARY(phist-olas STATIC ${PHIST_SRCS})

SET(TARGET_LL
  ${PHIST_LIBRARY}
  ghost
  hwloc
  phist_kernels_ghost
  phist_core
  phist_solvers
  phist_tools
  phist_precon
  #${LAPACK_LIBRARY}
  #${BLAS_LIBRARY}
  #${CFS_FORTRAN_LIBS}
  #matvec 
)

TARGET_LINK_LIBRARIES(phist-olas ${TARGET_LL})

