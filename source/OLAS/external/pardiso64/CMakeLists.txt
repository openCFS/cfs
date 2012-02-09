
IF(EXPLICIT_TEMPLATE_INSTANTIATION)
    SET(PARDISO64_SRCS pardiso64solver.cc)
ELSE(EXPLICIT_TEMPLATE_INSTANTIATION)
    SET(PARDISO64_SRCS Pardiso64Dummy.cc)
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(pardiso64-olas STATIC ${PARDISO64_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
TARGET_LINK_LIBRARIES(pardiso64-olas
  cfsgeneral
  datainout
  graph-olas
  logging
  matvec
  paramh
  solver-olas
  pardiso64-lib1
)

SET(DEPENDS
  "${LAPACK_LIBRARY}"
  "${BLAS_LIBRARY}"
  "${CFS_FORTRAN_LIBS}"
)

ADD_DEPENDENT_LIBRARY(pardiso64-lib "${PARDISO_LIBRARY}" "${DEPENDS}")
