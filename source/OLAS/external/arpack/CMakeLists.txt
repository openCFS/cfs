SET(ARPACK_SRCS
  ArpackEigenSolver.cc
  ArpackMatInterface.cc
  ArpackSolver.cc
  )

ADD_LIBRARY(arpack-olas STATIC ${ARPACK_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
SET(TARGET_LL
  ${ARPACK_LIBRARY}
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS}
  ${MUPARSER_LIBRARY}
  matvec
  )

ADD_DEPENDENCIES(arpack-olas arpack)

TARGET_LINK_LIBRARIES(arpack-olas ${TARGET_LL})


