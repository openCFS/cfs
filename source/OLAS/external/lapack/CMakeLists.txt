SET(LAPACK_SRCS
  LapackBaseMatrix.cc
  Lapack_LU.cc
  Lapack_LL.cc)

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(LAPACK_SRCS
    ${LAPACK_SRCS}
    LapackGBMatrix.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(lapack-olas STATIC ${LAPACK_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
SET(TARGET_LL
  matvec
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS})

TARGET_LINK_LIBRARIES(lapack-olas ${TARGET_LL})

