
INCLUDE_DIRECTORIES(${SUPERLU_INCLUDE_DIR})

SET(SUPERLU_SRCS SuperLUSolver.cc)

ADD_LIBRARY(superlu-olas STATIC ${SUPERLU_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
SET(TARGET_LL
  ${SUPERLU_LIBRARY}
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS})

TARGET_LINK_LIBRARIES(superlu-olas ${TARGET_LL})

ADD_DEPENDENCIES(superlu-olas cfsdeps)

