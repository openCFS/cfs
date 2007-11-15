SET(ILUPACK_SRCS Ilupack.cc)

ADD_LIBRARY(ilupack-olas STATIC ${ILUPACK_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
SET(TARGET_LL
  Dilupack
  Zilupack
  )

IF(BLAS_LIBRARY MATCHES "libmkl")
  SET(TARGET_LL
    ${TARGET_LL}
    blaslike_mkl
    )
ELSE(BLAS_LIBRARY MATCHES "libmkl")
  SET(TARGET_LL
    ${TARGET_LL}
    blaslike
    )
ENDIF(BLAS_LIBRARY MATCHES "libmkl")

SET(TARGET_LL
  ${TARGET_LL}
  sparspak
  amd
  ${PARDISO_LIBRARY}
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS})

TARGET_LINK_LIBRARIES(ilupack-olas ${TARGET_LL})

