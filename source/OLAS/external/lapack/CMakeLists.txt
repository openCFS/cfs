SET(LAPACK_SRCS
  lapacklu.cc
  lapackll.cc
  lapackbasematrix.cc
  olasf77mapping.cc
)

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(LAPACK_SRCS
    ${LAPACK_SRCS}
    lapackgbmatrix.cc
    )
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(lapack-olas STATIC ${LAPACK_SRCS})

# needed as solver-olas is cyclic dependent on us 
SET_PROPERTY(TARGET lapack-olas PROPERTY LINK_INTERFACE_MULTIPLICITY 3)

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
SET(TARGET_LL
  cfsgeneral
  graph-olas
  matvec
  paramh
  solver-olas
  utils
)

IF("${LAPACK_LIBRARY}" STREQUAL "")
ELSE("${LAPACK_LIBRARY}" STREQUAL "")

SET(DEPENDS
  "${BLAS_LIBRARY}"
  "${CFS_FORTRAN_LIBS}"
)

ADD_DEPENDENT_LIBRARY(lapack-lib "${LAPACK_LIBRARY}" "${DEPENDS}")

TARGET_LINK_LIBRARIES(lapack-olas ${TARGET_LL} lapack-lib1)

ENDIF("${LAPACK_LIBRARY}" STREQUAL "")

