
SET(UMFPACK_SRCS UMFPACKSolver.cc)

ADD_LIBRARY(umfpack-olas STATIC ${UMFPACK_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
SET(TARGET_LL
  ${UMFPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS})

TARGET_LINK_LIBRARIES(umfpack-olas ${TARGET_LL})

ADD_EXECUTABLE(umfpack_di_demo umfpack_di_demo.cc)

#IF(MINGW)
#  IF(CFS_BLAS_LAPACK STREQUAL MKL)
#    SET_TARGET_PROPERTIES(umfpack_di_demo PROPERTIES
#      LINK_FLAGS -Wl,--wrap,__chkstk
#      LINKER_LANGUAGE C)
#  ENDIF()
#ENDIF()

TARGET_LINK_LIBRARIES(umfpack_di_demo ${TARGET_LL})
