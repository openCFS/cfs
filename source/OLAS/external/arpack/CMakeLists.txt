SET(ARPACK_SRCS
  arpackEigensolver.cc
  arpackMatInterface.cc
  arpackSolver.cc
  )

ADD_LIBRARY(arpack-olas STATIC ${ARPACK_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
TARGET_LINK_LIBRARIES(arpack-olas
  cfsgeneral
  matvec
  paramh
  precond-olas
  solver-olas
  utils
  arpack-lib1
)

SET(DEPENDS
  "${LAPACK_LIBRARY}"
  "${BLAS_LIBRARY}"
  "${CFS_FORTRAN_LIBS}"
)

ADD_DEPENDENT_LIBRARY(arpack-lib "${ARPACK_LIBRARY}" "${DEPENDS}")

