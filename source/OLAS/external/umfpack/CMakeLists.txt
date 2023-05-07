
SET(UMFPACK_SRCS UMFPACKSolver.cc)

ADD_LIBRARY(umfpack-olas STATIC ${UMFPACK_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
set(TARGET_LL
  ${SUITESPARSE_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS})

# we had the case on centos6/gcc7 that clock_gettime@@GLIBC_2.2.5 is missing
if(UNIX AND NOT APPLE)
  set(TARGET_LL ${TARGET_LL} rt)
endif()

TARGET_LINK_LIBRARIES(umfpack-olas ${TARGET_LL})

ADD_DEPENDENCIES(umfpack-olas cfsdeps)
