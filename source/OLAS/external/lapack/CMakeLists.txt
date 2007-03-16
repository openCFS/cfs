SET(LAPACK_SRCS  olasf77mapping.cc lapackbasematrix.cc lapacklu.cc lapackll.cc)

ADD_LIBRARY(lapack-olas STATIC ${LAPACK_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
SET(TARGET_LL ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${CFS_FORTRAN_LIBS})

TARGET_LINK_LIBRARIES(lapack-olas ${TARGET_LL})

